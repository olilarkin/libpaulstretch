// Tests for the chunked offline render API (render_mono_chunked /
// render_stereo_chunked).
//
// These exist so the WebAssembly build can produce very long outputs without
// materialising the whole result in linear memory. renderMono returns the
// output in one buffer that must live in WASM memory twice over (the internal
// vector plus the returned copy), so a large render can exceed the wasm32 heap
// and abort. The chunked variants run the same DSP loop but deliver the result
// one bufsize() chunk at a time, keeping peak memory bounded.
//
// Part 1 checks behavioural parity with the whole-buffer path at a modest size.
// Part 2 checks the thing that actually matters — that a representative extreme
// render would blow the known wasm32 memory limits if buffered, while the
// chunked path's working set stays a tiny fraction of them. (The native host is
// 64-bit and can't reproduce the abort, so Part 2 reasons about the sizes via
// estimate_output_frames rather than allocating them.)
//
// We can't assert bit-identity against render_mono: each Stretcher instance
// bumps a static PRNG seed, so two renders get different phase randomization
// (same convention as streaming_test.cpp). Hence structural comparison.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "paulstretch/paulstretch.h"

namespace {

int failures = 0;

void check(bool cond, const std::string &label) {
	if (!cond) {
		++failures;
		std::cerr << "FAIL: " << label << '\n';
	} else {
		std::cout << "  ok: " << label << '\n';
	}
}

std::vector<float> make_sine(std::size_t frames, float sample_rate, float frequency) {
	std::vector<float> result(frames, 0.0f);
	const float tau = 6.28318530717958647692f;
	for (std::size_t i = 0; i < frames; i++) {
		float phase = tau * frequency * static_cast<float>(i) / sample_rate;
		result[i] = 0.5f * std::sin(phase);
	}
	return result;
}

struct Stats {
	std::size_t frames = 0;
	float peak = 0.0f;
	double rms = 0.0;
	bool all_finite = true;
};

Stats analyze(const std::vector<float> &v) {
	Stats s;
	s.frames = v.size();
	double sumsq = 0.0;
	for (float x : v) {
		if (!std::isfinite(x)) { s.all_finite = false; continue; }
		float a = std::fabs(x);
		if (a > s.peak) s.peak = a;
		sumsq += static_cast<double>(x) * static_cast<double>(x);
	}
	s.rms = v.empty() ? 0.0 : std::sqrt(sumsq / static_cast<double>(v.size()));
	return s;
}

// Same contract streaming_test uses: equal length, finite, peak within 2×, RMS
// within 15% (phase randomness differs but averages out).
void check_structural(const std::vector<float> &a, const std::vector<float> &b,
                      const std::string &label) {
	auto sa = analyze(a);
	auto sb = analyze(b);
	check(sa.all_finite, label + ": chunked finite");
	check(sb.all_finite, label + ": reference finite");
	check(a.size() == b.size(), label + ": same length");
	if (sa.peak > 0 && sb.peak > 0) {
		const double ratio = sa.peak / sb.peak;
		check(ratio > 0.5 && ratio < 2.0, label + ": peaks within 2×");
	}
	if (sa.rms > 0 && sb.rms > 0) {
		const double rel = std::fabs(sa.rms - sb.rms) / std::max(sa.rms, sb.rms);
		check(rel < 0.15, label + ": RMS within 15%");
	}
}

double gib(std::uint64_t bytes) { return static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0); }

// ── Part 1: behavioural parity at a modest, CI-fast size ───────────────────
void test_parity() {
	std::cout << "\n[parity] chunked vs whole-buffer render\n";

	const paulstretch::RenderOptions options{
		.stretch = 32.0f,
		.fft_size = 4096,
		.sample_rate = 44100.0f,
		.window = paulstretch::Window::Hann,
		.onset_detection_sensitivity = 0.0f,
	};

	paulstretch::OfflineRenderer renderer(options);
	const int bufsize = static_cast<int>(renderer.options().fft_size);

	const auto left_in = make_sine(44100, options.sample_rate, 220.0f);  // 1 s
	const auto right_in = make_sine(44100, options.sample_rate, 277.0f);

	// Mono.
	const auto mono_ref = renderer.render_mono(left_in);
	std::vector<float> mono_chunked;
	mono_chunked.reserve(mono_ref.size());
	std::size_t mono_chunks = 0;
	bool mono_uniform = true;
	renderer.render_mono_chunked(left_in, [&](const float *data, int frames) {
		mono_chunks++;
		if (frames != bufsize) mono_uniform = false;
		mono_chunked.insert(mono_chunked.end(), data, data + frames);
	});
	std::cout << "  mono: " << mono_ref.size() << " frames in " << mono_chunks << " chunks\n";
	check(!mono_ref.empty(), "render_mono produced output");
	check(mono_uniform, "every mono chunk is exactly bufsize frames");
	check(mono_chunks > 100, "many mono chunks (per-chunk memory is bounded)");
	check_structural(mono_chunked, mono_ref, "mono");

	// Stereo.
	const auto stereo_ref = renderer.render_stereo(left_in, right_in);
	std::vector<float> stereo_l, stereo_r;
	stereo_l.reserve(stereo_ref.left.size());
	stereo_r.reserve(stereo_ref.right.size());
	bool stereo_uniform = true;
	renderer.render_stereo_chunked(left_in, right_in, [&](const float *l, const float *r, int frames) {
		if (frames != bufsize) stereo_uniform = false;
		stereo_l.insert(stereo_l.end(), l, l + frames);
		stereo_r.insert(stereo_r.end(), r, r + frames);
	});
	check(!stereo_ref.left.empty(), "render_stereo produced output");
	check(stereo_uniform, "every stereo chunk is exactly bufsize frames");
	check_structural(stereo_l, stereo_ref.left, "stereo L");
	check_structural(stereo_r, stereo_ref.right, "stereo R");

	// Empty input is a no-op (sink never fires).
	std::size_t empty_calls = 0;
	renderer.render_mono_chunked({}, [&](const float *, int) { empty_calls++; });
	check(empty_calls == 0, "chunked render of empty input does not call the sink");
}

// ── Part 2: the chunked path stays under the wasm32 memory limits ──────────
void test_memory_bounds() {
	std::cout << "\n[limits] buffered footprint vs wasm32 heap caps\n";

	// wasm32 linear-memory ceilings the WASM build runs against.
	constexpr std::uint64_t kFloat = sizeof(float);
	constexpr std::uint64_t kWasmDefaultCap = 2ull * 1024 * 1024 * 1024; // Emscripten default MAXIMUM_MEMORY
	constexpr std::uint64_t kWasmMaxHeap = 4ull * 1024 * 1024 * 1024;    // wasm32 hard maximum

	// A representative extreme job: a 15-second clip stretched 500× (~2 hours of
	// output). Sizes come from estimate_output_frames — nothing is allocated.
	const paulstretch::RenderOptions options{
		.stretch = 500.0f,
		.fft_size = 4096,
		.sample_rate = 44100.0f,
		.window = paulstretch::Window::Hann,
		.onset_detection_sensitivity = 0.0f,
	};
	paulstretch::OfflineRenderer renderer(options);

	const std::size_t input_frames = static_cast<std::size_t>(15 * 44100);
	const std::uint64_t out_frames = renderer.estimate_output_frames(input_frames);
	const std::uint64_t per_channel = out_frames * kFloat;

	// What renderMono / renderStereo actually need live at once:
	//   mono   — internal vector + the returned JS copy            = 2× a channel
	//   stereo — both channel vectors, plus one channel's JS copy  = 3× a channel
	const std::uint64_t mono_peak = per_channel * 2;
	const std::uint64_t stereo_peak = per_channel * 3;

	// What the chunked path needs live at once: a single bufsize chunk.
	const std::uint64_t chunked_working_set =
		static_cast<std::uint64_t>(renderer.options().fft_size) * kFloat;

	std::cout << "  output ~" << out_frames << " frames (" << gib(per_channel) << " GiB/channel)\n";
	std::cout << "  buffered peak: mono " << gib(mono_peak) << " GiB, stereo "
	          << gib(stereo_peak) << " GiB\n";
	std::cout << "  chunked working set: " << chunked_working_set << " bytes\n";

	// Buffered rendering of this job would have aborted under Emscripten's old
	// 2 GiB default heap cap — both mono and stereo overflow it.
	check(mono_peak > kWasmDefaultCap,
	      "buffered mono render exceeds the 2 GiB default heap cap");
	check(stereo_peak > kWasmDefaultCap,
	      "buffered stereo render exceeds the 2 GiB default heap cap");
	// Raising the cap to the wasm32 maximum (4 GiB) is what lets this particular
	// job through on the buffered path.
	check(stereo_peak < kWasmMaxHeap,
	      "...but fits within the raised 4 GiB cap");

	// The chunked path's footprint is a negligible fraction of the cap, so it
	// completes regardless of how long the output is.
	check(chunked_working_set < (1ull << 20),
	      "chunked working set is under 1 MiB");
	check(chunked_working_set * 1000 < kWasmDefaultCap,
	      "chunked working set is <0.1% of the heap cap");

	// A bigger job (30 s × 500×) overflows even the 4 GiB maximum on the *mono*
	// buffered path — so raising the cap is not enough, and only the chunked
	// path (which is independent of output length) stays under the limit.
	paulstretch::OfflineRenderer big({
		.stretch = 500.0f,
		.fft_size = 4096,
		.sample_rate = 48000.0f,
		.window = paulstretch::Window::Hann,
		.onset_detection_sensitivity = 0.0f,
	});
	const std::uint64_t big_per_channel =
		static_cast<std::uint64_t>(big.estimate_output_frames(30 * 48000)) * kFloat;
	std::cout << "  bigger job buffered mono peak: " << gib(big_per_channel * 2) << " GiB\n";
	check(big_per_channel * 2 > kWasmMaxHeap,
	      "a bigger job exceeds even the 4 GiB wasm32 maximum (chunking required)");
}

} // namespace

int main() {
	std::cout << "backend=" << paulstretch::fft_backend_name()
	          << " simd=" << paulstretch::fft_simd_arch()
	          << " width=" << paulstretch::fft_simd_size() << '\n';

	test_parity();
	test_memory_bounds();

	if (failures == 0) {
		std::cout << "\nOK\n";
		return EXIT_SUCCESS;
	}
	std::cerr << "\n" << failures << " check(s) failed\n";
	return EXIT_FAILURE;
}
