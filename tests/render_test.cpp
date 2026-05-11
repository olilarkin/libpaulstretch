// Render-properties test: drives the OfflineRenderer through a sweep of
// configurations against plenty_of_unknown.wav and asserts numerical
// invariants. Path to the WAV is passed as argv[1] (set by CMake).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

#include "paulstretch/paulstretch.h"

namespace {

struct Wav {
    std::vector<float> left;
    std::vector<float> right;
    unsigned int channels = 0;
    unsigned int sample_rate = 0;
};

bool load_wav(const std::string &path, Wav &out) {
    drwav wav;
    if (!drwav_init_file(&wav, path.c_str(), nullptr)) return false;
    out.channels = wav.channels;
    out.sample_rate = wav.sampleRate;
    const std::size_t frames = static_cast<std::size_t>(wav.totalPCMFrameCount);
    std::vector<float> interleaved(frames * wav.channels);
    drwav_read_pcm_frames_f32(&wav, frames, interleaved.data());
    drwav_uninit(&wav);
    out.left.resize(frames);
    out.right.resize(frames);
    if (wav.channels == 1) {
        out.left = std::move(interleaved);
        out.right = out.left;
    } else {
        for (std::size_t i = 0; i < frames; i++) {
            out.left[i] = interleaved[i * 2];
            out.right[i] = interleaved[i * 2 + 1];
        }
    }
    return true;
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

int failures = 0;

void check(bool cond, const std::string &label) {
    if (!cond) {
        ++failures;
        std::cerr << "FAIL: " << label << '\n';
    } else {
        std::cout << "  ok: " << label << '\n';
    }
}

// FNV-1a hash of the raw float bit-pattern — cheap and stable across runs.
std::uint64_t hash_floats(const std::vector<float> &v) {
    std::uint64_t h = 1469598103934665603ull;
    for (float x : v) {
        std::uint32_t bits;
        std::memcpy(&bits, &x, sizeof bits);
        for (int b = 0; b < 4; b++) {
            h ^= (bits >> (b * 8)) & 0xff;
            h *= 1099511628211ull;
        }
    }
    return h;
}

struct Case {
    const char *name;
    paulstretch::RenderOptions opts;
};

void run_mono_case(const Case &c, const Wav &input) {
    std::cout << "\n[case] " << c.name << "\n";
    paulstretch::OfflineRenderer renderer(c.opts);
    const auto &sanitized = renderer.options();
    std::cout << "  fft_size=" << sanitized.fft_size
              << " stretch=" << sanitized.stretch
              << " onset=" << sanitized.onset_detection_sensitivity
              << '\n';

    const std::size_t estimate = renderer.estimate_output_frames(input.left.size());
    auto out = renderer.render_mono(input.left);
    auto s = analyze(out);

    std::cout << "  frames=" << s.frames
              << " estimate=" << estimate
              << " peak=" << s.peak
              << " rms=" << s.rms
              << " hash=" << std::hex << hash_floats(out) << std::dec
              << '\n';

    check(s.all_finite, std::string(c.name) + ": all samples finite");
    check(s.frames > input.left.size(),
          std::string(c.name) + ": stretched output longer than input");
    // estimate is an upper bound used for buffer sizing; actual must not exceed it.
    check(s.frames <= estimate,
          std::string(c.name) + ": frames <= estimate_output_frames");
    check(s.peak <= 4.0f,
          std::string(c.name) + ": peak within reasonable bound");
    check(s.peak > 1e-4f,
          std::string(c.name) + ": output is not silent");
}

void run_stereo_case(const Case &c, const Wav &input) {
    std::cout << "\n[stereo case] " << c.name << "\n";
    paulstretch::OfflineRenderer renderer(c.opts);
    auto out = renderer.render_stereo(input.left, input.right);

    auto sL = analyze(out.left);
    auto sR = analyze(out.right);

    std::cout << "  L: frames=" << sL.frames << " peak=" << sL.peak
              << " rms=" << sL.rms << " hash=" << std::hex
              << hash_floats(out.left) << std::dec << '\n';
    std::cout << "  R: frames=" << sR.frames << " peak=" << sR.peak
              << " rms=" << sR.rms << " hash=" << std::hex
              << hash_floats(out.right) << std::dec << '\n';

    check(sL.all_finite && sR.all_finite,
          std::string(c.name) + ": stereo all finite");
    check(sL.frames == sR.frames,
          std::string(c.name) + ": stereo channels same length");
    // Per-channel phase randomization should make the channels differ for
    // a non-degenerate stereo input (different seeds per Stretcher).
    check(hash_floats(out.left) != hash_floats(out.right),
          std::string(c.name) + ": stereo channels independent");
    // But peaks/RMS should be in the same ballpark since the source content
    // is similar between channels.
    check(std::fabs(sL.rms - sR.rms) / std::max(sL.rms, sR.rms) < 0.5,
          std::string(c.name) + ": stereo channel RMS within 50%");
}

void run_envelope_case(const Wav &input) {
    std::cout << "\n[envelope case] linear ramp 2x -> 8x\n";
    paulstretch::RenderOptions opts;
    opts.stretch = 1.0f; // overridden by envelope
    opts.fft_size = 2048;
    opts.sample_rate = static_cast<float>(input.sample_rate);
    opts.window = paulstretch::Window::Hann;

    paulstretch::OfflineRenderer renderer_flat(opts);
    paulstretch::OfflineRenderer renderer_env(opts);
    renderer_env.set_stretch_envelope({{0.0f, 2.0f}, {1.0f, 8.0f}});

    auto flat = renderer_flat.render_mono(input.left);
    auto env  = renderer_env.render_mono(input.left);

    const auto flat_stats = analyze(flat);
    const auto env_stats  = analyze(env);

    std::cout << "  flat frames=" << flat_stats.frames
              << " env frames=" << env_stats.frames << '\n';

    check(env_stats.all_finite, "envelope: output finite");
    // Ramp 2 -> 8 should average roughly 5x; flat is 1x. Output must be longer.
    check(env_stats.frames > flat_stats.frames * 3,
          "envelope: ramp output longer than flat 1x by >3x");
    check(env_stats.frames < flat_stats.frames * 10,
          "envelope: ramp output not absurdly long");
}

void run_window_sweep(const Wav &input) {
    std::cout << "\n[window sweep]\n";
    const paulstretch::Window windows[] = {
        paulstretch::Window::Rectangular,
        paulstretch::Window::Hamming,
        paulstretch::Window::Hann,
        paulstretch::Window::Blackman,
        paulstretch::Window::BlackmanHarris,
    };
    const char *names[] = {
        "rectangular", "hamming", "hann", "blackman", "blackman-harris",
    };
    std::vector<std::uint64_t> hashes;
    for (int i = 0; i < 5; i++) {
        paulstretch::RenderOptions opts;
        opts.stretch = 4.0f;
        opts.fft_size = 2048;
        opts.sample_rate = static_cast<float>(input.sample_rate);
        opts.window = windows[i];
        paulstretch::OfflineRenderer r(opts);
        auto out = r.render_mono(input.left);
        auto s = analyze(out);
        auto h = hash_floats(out);
        hashes.push_back(h);
        std::cout << "  " << names[i] << ": frames=" << s.frames
                  << " peak=" << s.peak << " hash=" << std::hex << h
                  << std::dec << '\n';
        check(s.all_finite, std::string("window ") + names[i] + ": finite");
    }
    // Each window should give a distinct output.
    for (std::size_t i = 0; i < hashes.size(); i++) {
        for (std::size_t j = i + 1; j < hashes.size(); j++) {
            check(hashes[i] != hashes[j],
                  std::string("windows differ: ") + names[i] + " vs " + names[j]);
        }
    }
}

void run_determinism(const Wav &input) {
    std::cout << "\n[determinism]\n";
    paulstretch::RenderOptions opts;
    opts.stretch = 4.0f;
    opts.fft_size = 2048;
    opts.sample_rate = static_cast<float>(input.sample_rate);
    opts.window = paulstretch::Window::Hann;
    // Two renderers in the SAME process share start_rand_seed_, so seed
    // advances between them — outputs should differ. (This is documented
    // behavior of the algorithm's per-Stretcher seed counter.)
    paulstretch::OfflineRenderer r1(opts);
    paulstretch::OfflineRenderer r2(opts);
    auto a = r1.render_mono(input.left);
    auto b = r2.render_mono(input.left);
    check(hash_floats(a) != hash_floats(b),
          "different renderers in one process produce different output (per-instance seed)");
    check(a.size() == b.size(),
          "different renderers in one process produce same-length output");
}

} // namespace

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <input.wav>\n";
        return EXIT_FAILURE;
    }
    Wav input;
    if (!load_wav(argv[1], input)) {
        std::cerr << "error: cannot open " << argv[1] << '\n';
        return EXIT_FAILURE;
    }
    std::cout << "loaded " << input.left.size() << " frames, "
              << input.channels << "ch @ " << input.sample_rate << " Hz\n";
    std::cout << "backend=" << paulstretch::fft_backend_name()
              << " simd=" << paulstretch::fft_simd_arch()
              << " width=" << paulstretch::fft_simd_size() << '\n';

    const float sr = static_cast<float>(input.sample_rate);

    const Case mono_cases[] = {
        {"stretch=2 fft=1024 hann",
         {.stretch = 2.0f,  .fft_size = 1024, .sample_rate = sr,
          .window = paulstretch::Window::Hann}},
        {"stretch=8 fft=4096 hann",
         {.stretch = 8.0f,  .fft_size = 4096, .sample_rate = sr,
          .window = paulstretch::Window::Hann}},
        {"stretch=16 fft=8192 blackman-harris",
         {.stretch = 16.0f, .fft_size = 8192, .sample_rate = sr,
          .window = paulstretch::Window::BlackmanHarris}},
        {"stretch=4 fft=2048 hann onset=0.5",
         {.stretch = 4.0f,  .fft_size = 2048, .sample_rate = sr,
          .window = paulstretch::Window::Hann,
          .onset_detection_sensitivity = 0.5f}},
    };

    for (const auto &c : mono_cases) run_mono_case(c, input);

    const Case stereo_cases[] = {
        {"stereo stretch=4 fft=2048 hann",
         {.stretch = 4.0f, .fft_size = 2048, .sample_rate = sr,
          .window = paulstretch::Window::Hann}},
        {"stereo stretch=12 fft=8192 blackman",
         {.stretch = 12.0f, .fft_size = 8192, .sample_rate = sr,
          .window = paulstretch::Window::Blackman}},
    };
    for (const auto &c : stereo_cases) run_stereo_case(c, input);

    run_window_sweep(input);
    run_envelope_case(input);
    run_determinism(input);

    std::cout << "\n" << (failures ? "FAILED" : "PASSED")
              << ": " << failures << " failure(s)\n";
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
