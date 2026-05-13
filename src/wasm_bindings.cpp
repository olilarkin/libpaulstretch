#include <algorithm>
#include <stdexcept>
#include <vector>

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "paulstretch/paulstretch.h"

namespace {

std::vector<float> from_js_array(const emscripten::val &input) {
	const std::size_t length = input["length"].as<std::size_t>();
	std::vector<float> samples(length, 0.0f);
	emscripten::val view = emscripten::val(emscripten::typed_memory_view(length, samples.data()));
	view.call<void>("set", input);
	return samples;
}

emscripten::val to_js_float32_array(const std::vector<float> &input) {
	emscripten::val output = emscripten::val::global("Float32Array").new_(input.size());
	output.call<void>("set", emscripten::val(emscripten::typed_memory_view(input.size(), input.data())));
	return output;
}

template <typename T>
T get_or(const emscripten::val &obj, const char *key, T fallback) {
	emscripten::val v = obj[key];
	if (v.isNull() || v.isUndefined()) return fallback;
	return v.as<T>();
}

paulstretch::ProcessOptions process_options_from_js(const emscripten::val &obj) {
	paulstretch::ProcessOptions options;
	options.pitch_shift_enabled = get_or(obj, "pitchShiftEnabled", options.pitch_shift_enabled);
	options.pitch_shift_cents = get_or(obj, "pitchShiftCents", options.pitch_shift_cents);

	options.octave_enabled = get_or(obj, "octaveEnabled", options.octave_enabled);
	options.octave_minus2 = get_or(obj, "octaveMinus2", options.octave_minus2);
	options.octave_minus1 = get_or(obj, "octaveMinus1", options.octave_minus1);
	options.octave_0 = get_or(obj, "octave0", options.octave_0);
	options.octave_plus1 = get_or(obj, "octavePlus1", options.octave_plus1);
	options.octave_plus15 = get_or(obj, "octavePlus15", options.octave_plus15);
	options.octave_plus2 = get_or(obj, "octavePlus2", options.octave_plus2);

	options.frequency_shift_enabled = get_or(obj, "frequencyShiftEnabled", options.frequency_shift_enabled);
	options.frequency_shift_hz = get_or(obj, "frequencyShiftHz", options.frequency_shift_hz);

	options.compressor_enabled = get_or(obj, "compressorEnabled", options.compressor_enabled);
	options.compressor_power = get_or(obj, "compressorPower", options.compressor_power);

	options.filter_enabled = get_or(obj, "filterEnabled", options.filter_enabled);
	options.filter_low_hz = get_or(obj, "filterLowHz", options.filter_low_hz);
	options.filter_high_hz = get_or(obj, "filterHighHz", options.filter_high_hz);
	options.filter_high_damp = get_or(obj, "filterHighDamp", options.filter_high_damp);
	options.filter_stop = get_or(obj, "filterStop", options.filter_stop);

	options.harmonics_enabled = get_or(obj, "harmonicsEnabled", options.harmonics_enabled);
	options.harmonics_frequency_hz = get_or(obj, "harmonicsFrequencyHz", options.harmonics_frequency_hz);
	options.harmonics_bandwidth_cents = get_or(obj, "harmonicsBandwidthCents", options.harmonics_bandwidth_cents);
	options.harmonics_count = get_or(obj, "harmonicsCount", options.harmonics_count);
	options.harmonics_gauss = get_or(obj, "harmonicsGauss", options.harmonics_gauss);

	options.spread_enabled = get_or(obj, "spreadEnabled", options.spread_enabled);
	options.spread_bandwidth = get_or(obj, "spreadBandwidth", options.spread_bandwidth);

	options.tonal_noise_enabled = get_or(obj, "tonalNoiseEnabled", options.tonal_noise_enabled);
	options.tonal_noise_preserve = get_or(obj, "tonalNoisePreserve", options.tonal_noise_preserve);
	options.tonal_noise_bandwidth = get_or(obj, "tonalNoiseBandwidth", options.tonal_noise_bandwidth);

	options.arbitrary_filter_enabled = get_or(obj, "arbitraryFilterEnabled", options.arbitrary_filter_enabled);
	return options;
}

std::vector<paulstretch::Breakpoint> breakpoints_from_js(
	const emscripten::val &xs,
	const emscripten::val &ys) {
	const int n = std::min(xs["length"].as<int>(), ys["length"].as<int>());
	std::vector<paulstretch::Breakpoint> points(n);
	for (int i = 0; i < n; i++) {
		points[i].position = xs[i].as<float>();
		points[i].value = ys[i].as<float>();
	}
	return points;
}

class WasmOfflineRenderer {
public:
	WasmOfflineRenderer()
		: renderer_() {}

	WasmOfflineRenderer(
		float stretch,
		int fft_size,
		float sample_rate,
		paulstretch::Window window,
		float onset_detection_sensitivity)
		: renderer_({
			stretch,
			fft_size,
			sample_rate,
			window,
			onset_detection_sensitivity,
		}) {}

	emscripten::val renderMono(const emscripten::val &input) const {
		return to_js_float32_array(renderer_.render_mono(from_js_array(input)));
	}

	emscripten::val renderStereo(const emscripten::val &left, const emscripten::val &right) const {
		const auto stereo = renderer_.render_stereo(from_js_array(left), from_js_array(right));
		emscripten::val result = emscripten::val::object();
		result.set("left", to_js_float32_array(stereo.left));
		result.set("right", to_js_float32_array(stereo.right));
		return result;
	}

	void setStretchEnvelope(const emscripten::val &xs, const emscripten::val &ys) {
		const int n = std::min(xs["length"].as<int>(), ys["length"].as<int>());
		std::vector<paulstretch::Breakpoint> envelope(n);
		for (int i = 0; i < n; i++) {
			envelope[i].position = xs[i].as<float>();
			envelope[i].value = ys[i].as<float>();
		}
		renderer_.set_stretch_envelope(std::move(envelope));
	}

	void clearStretchEnvelope() {
		renderer_.clear_stretch_envelope();
	}

	void setProcessOptions(const emscripten::val &options) {
		renderer_.set_process_options(process_options_from_js(options));
	}

	void setArbitraryFilter(const emscripten::val &xs, const emscripten::val &ys) {
		renderer_.set_arbitrary_filter(breakpoints_from_js(xs, ys));
	}

	void clearArbitraryFilter() {
		renderer_.clear_arbitrary_filter();
	}

	std::size_t estimateOutputFrames(std::size_t input_frames) const {
		return renderer_.estimate_output_frames(input_frames);
	}

private:
	paulstretch::OfflineRenderer renderer_;
};

// Wasm-facing wrapper around StreamingStretcher. Each step() call allocates a
// new Float32Array for the output chunk — for typical realtime use that's one
// allocation per ~bufsize() output frames (~10/sec at bufsize=4096), well below
// any concerning GC pressure. A future zero-copy variant could expose a stable
// typed_memory_view into a fixed internal buffer, at the cost of invalidating
// the view whenever the wasm heap grows (ALLOW_MEMORY_GROWTH=1).
class WasmStreamingStretcher {
public:
	WasmStreamingStretcher(
		float stretch,
		int fft_size,
		float sample_rate,
		paulstretch::Window window,
		float onset_detection_sensitivity)
		: inner_({
			stretch,
			fft_size,
			sample_rate,
			window,
			onset_detection_sensitivity,
		}),
		  out_buf_(inner_.bufsize(), 0.0f) {}

	int bufsize() const { return inner_.bufsize(); }
	int maxInputChunk() const { return inner_.max_input_chunk(); }
	int nextInputSize() const { return inner_.next_input_size(); }
	int skipAfterStep() const { return inner_.skip_after_step(); }

	// Push one chunk of input (or null/undefined for the no-input case),
	// produce one chunk of output. `position_pct` is the input cursor as a
	// percent 0..100 for envelope evaluation. Returns an object:
	//   { output: Float32Array(bufsize()), onset: number }
	emscripten::val step(const emscripten::val &input, float position_pct) {
		return stepImpl(input, position_pct, true);
	}

	emscripten::val stepWithoutOnsetFeedback(const emscripten::val &input, float position_pct) {
		return stepImpl(input, position_pct, false);
	}

	void applyOnset(float onset) { inner_.apply_onset(onset); }

private:
	emscripten::val stepImpl(const emscripten::val &input, float position_pct, bool apply_onset) {
		std::vector<float> in_local;
		const float *in_ptr = nullptr;
		if (!input.isNull() && !input.isUndefined()) {
			in_local = from_js_array(input);
			if (!in_local.empty()) in_ptr = in_local.data();
		}
		const float onset = apply_onset
			? inner_.step(in_ptr, position_pct, out_buf_.data())
			: inner_.step_without_onset_feedback(in_ptr, position_pct, out_buf_.data());

		emscripten::val result = emscripten::val::object();
		result.set("output", to_js_float32_array(out_buf_));
		result.set("onset", onset);
		return result;
	}

public:
	void setStretchEnvelope(const emscripten::val &xs, const emscripten::val &ys) {
		inner_.set_stretch_envelope(breakpoints_from_js(xs, ys));
	}

	void clearStretchEnvelope() { inner_.clear_stretch_envelope(); }

	void setProcessOptions(const emscripten::val &options) {
		inner_.set_process_options(process_options_from_js(options));
	}

	void setArbitraryFilter(const emscripten::val &xs, const emscripten::val &ys) {
		inner_.set_arbitrary_filter(breakpoints_from_js(xs, ys));
	}

	void clearArbitraryFilter() { inner_.clear_arbitrary_filter(); }

	void setStretchFactor(float stretch) { inner_.set_stretch_factor(stretch); }

	void setOnsetDetectionSensitivity(float s) {
		inner_.set_onset_detection_sensitivity(s);
	}

	void reset() { inner_.reset(); }

private:
	paulstretch::StreamingStretcher inner_;
	std::vector<float> out_buf_;
};

class WasmBinauralBeatsProcessor {
public:
	explicit WasmBinauralBeatsProcessor(float sample_rate)
		: inner_(sample_rate) {}

	void setOptions(const emscripten::val &obj) {
		paulstretch::BinauralBeatsOptions options;
		options.enabled = get_or(obj, "enabled", options.enabled);
		options.stereo_mode = static_cast<paulstretch::BinauralStereoMode>(
			get_or(obj, "stereoMode", static_cast<int>(options.stereo_mode)));
		options.mono = get_or(obj, "mono", options.mono);
		options.beat_frequency_hz = get_or(obj, "beatFrequencyHz", options.beat_frequency_hz);
		inner_.set_options(options);
	}

	void setFrequencyEnvelope(const emscripten::val &xs, const emscripten::val &ys) {
		inner_.set_frequency_envelope(breakpoints_from_js(xs, ys));
	}

	void clearFrequencyEnvelope() {
		inner_.clear_frequency_envelope();
	}

	emscripten::val process(const emscripten::val &left, const emscripten::val &right, float position_pct) {
		std::vector<float> l = from_js_array(left);
		std::vector<float> r = from_js_array(right);
		if (l.size() != r.size()) throw std::invalid_argument("left and right channel lengths must match");
		inner_.process(l.data(), r.data(), static_cast<int>(l.size()), position_pct);

		emscripten::val result = emscripten::val::object();
		result.set("left", to_js_float32_array(l));
		result.set("right", to_js_float32_array(r));
		return result;
	}

	void reset() { inner_.reset(); }

private:
	paulstretch::BinauralBeatsProcessor inner_;
};

} // namespace

EMSCRIPTEN_BINDINGS(paulstretch) {
	emscripten::enum_<paulstretch::Window>("Window")
		.value("Rectangular", paulstretch::Window::Rectangular)
		.value("Hamming", paulstretch::Window::Hamming)
		.value("Hann", paulstretch::Window::Hann)
		.value("Blackman", paulstretch::Window::Blackman)
		.value("BlackmanHarris", paulstretch::Window::BlackmanHarris);

	emscripten::enum_<paulstretch::BinauralStereoMode>("BinauralStereoMode")
		.value("LeftRight", paulstretch::BinauralStereoMode::LeftRight)
		.value("RightLeft", paulstretch::BinauralStereoMode::RightLeft)
		.value("Symmetric", paulstretch::BinauralStereoMode::Symmetric);

	emscripten::class_<WasmOfflineRenderer>("OfflineRenderer")
		.constructor<>()
		.constructor<float, int, float, paulstretch::Window, float>()
		.function("estimateOutputFrames", &WasmOfflineRenderer::estimateOutputFrames)
		.function("renderMono", &WasmOfflineRenderer::renderMono)
		.function("renderStereo", &WasmOfflineRenderer::renderStereo)
		.function("setStretchEnvelope", &WasmOfflineRenderer::setStretchEnvelope)
		.function("clearStretchEnvelope", &WasmOfflineRenderer::clearStretchEnvelope)
		.function("setProcessOptions", &WasmOfflineRenderer::setProcessOptions)
		.function("setArbitraryFilter", &WasmOfflineRenderer::setArbitraryFilter)
		.function("clearArbitraryFilter", &WasmOfflineRenderer::clearArbitraryFilter);

	emscripten::class_<WasmStreamingStretcher>("StreamingStretcher")
		.constructor<float, int, float, paulstretch::Window, float>()
		.function("bufsize", &WasmStreamingStretcher::bufsize)
		.function("maxInputChunk", &WasmStreamingStretcher::maxInputChunk)
		.function("nextInputSize", &WasmStreamingStretcher::nextInputSize)
		.function("skipAfterStep", &WasmStreamingStretcher::skipAfterStep)
		.function("step", &WasmStreamingStretcher::step)
		.function("stepWithoutOnsetFeedback",
		          &WasmStreamingStretcher::stepWithoutOnsetFeedback)
		.function("applyOnset", &WasmStreamingStretcher::applyOnset)
		.function("setStretchEnvelope", &WasmStreamingStretcher::setStretchEnvelope)
		.function("clearStretchEnvelope", &WasmStreamingStretcher::clearStretchEnvelope)
		.function("setProcessOptions", &WasmStreamingStretcher::setProcessOptions)
		.function("setArbitraryFilter", &WasmStreamingStretcher::setArbitraryFilter)
		.function("clearArbitraryFilter", &WasmStreamingStretcher::clearArbitraryFilter)
		.function("setStretchFactor", &WasmStreamingStretcher::setStretchFactor)
		.function("setOnsetDetectionSensitivity",
		          &WasmStreamingStretcher::setOnsetDetectionSensitivity)
		.function("reset", &WasmStreamingStretcher::reset);

	emscripten::class_<WasmBinauralBeatsProcessor>("BinauralBeatsProcessor")
		.constructor<float>()
		.function("setOptions", &WasmBinauralBeatsProcessor::setOptions)
		.function("setFrequencyEnvelope", &WasmBinauralBeatsProcessor::setFrequencyEnvelope)
		.function("clearFrequencyEnvelope", &WasmBinauralBeatsProcessor::clearFrequencyEnvelope)
		.function("process", &WasmBinauralBeatsProcessor::process)
		.function("reset", &WasmBinauralBeatsProcessor::reset);

	emscripten::function("fftBackendName", &paulstretch::fft_backend_name);
	emscripten::function("fftSimdArch", &paulstretch::fft_simd_arch);
	emscripten::function("fftSimdSize", &paulstretch::fft_simd_size);
}
