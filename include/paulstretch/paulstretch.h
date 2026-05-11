#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace paulstretch {

enum class Window { Rectangular, Hamming, Hann, Blackman, BlackmanHarris };

struct Breakpoint {
    float position; // normalized time position in the input (0-1)
    float value;    // stretch multiplier at this position
};

struct RenderOptions {
    float stretch = 8.0f;
    int fft_size = 4096;
    float sample_rate = 44100.0f;
    Window window = Window::Hann;
    float onset_detection_sensitivity = 0.0f;
};

struct StereoBuffer {
    std::vector<float> left;
    std::vector<float> right;
};

std::string fft_backend_name();
std::string fft_simd_arch();
int fft_simd_size();

// ── StreamingStretcher ─────────────────────────────────────────────────────
//
// Block-based, push/pull stretching primitive. Designed for realtime use:
// the host (e.g. an AudioWorklet or Web Worker) feeds input in chunks and
// drains one output chunk of `bufsize()` frames per step, with no full-output
// buffering inside the library. The protocol is:
//
//   * Construct with sanitized options.
//   * Call next_input_size() to learn how many input frames step() wants.
//     On the very first call this is max_input_chunk() (initial fill);
//     thereafter it alternates between 0 and bufsize() depending on stretch
//     factor and onset detection.
//   * If non-zero, gather that many input frames from your source (zero-pad
//     if the source ran out). Pass them to step() along with the current
//     input cursor position as a percent 0..100 (used to evaluate the
//     stretch envelope).
//   * step() writes bufsize() output frames to the supplied output buffer
//     and returns the detected onset value.
//   * Advance your input cursor by skip_after_step() frames after step()
//     (in addition to whatever you consumed). For low stretches with onset
//     detection this can be substantial.
//
// The OfflineRenderer below is implemented on top of StreamingStretcher.
class StreamingStretcher {
public:
    explicit StreamingStretcher(RenderOptions options);
    ~StreamingStretcher();
    StreamingStretcher(const StreamingStretcher &) = delete;
    StreamingStretcher &operator=(const StreamingStretcher &) = delete;

    // Sanitized options actually in use (fft_size may have been rounded to a
    // backend-compatible size). Read this if you need the effective fft_size.
    const RenderOptions &options() const;

    // Frames of input step() wants on the next call: max_input_chunk() on the
    // first call (initial fill), then 0 or bufsize() thereafter.
    int next_input_size() const;

    // Output chunk size (sanitized fft_size). Each step() writes this many
    // frames into the caller's output buffer.
    int bufsize() const;

    // Largest single input chunk: 3 * bufsize(). The first step() call
    // expects this many frames.
    int max_input_chunk() const;

    // After step() returns, this is the number of additional input frames
    // the caller must skip before the next step() (in addition to the
    // next_input_size() it will request). Non-zero with onset detection
    // and at low stretch factors.
    int skip_after_step() const;

    // Run one step. `input` points to next_input_size() frames (may be
    // nullptr if next_input_size()==0). `position_pct` is the input
    // cursor as a percent 0..100 used to evaluate the stretch envelope.
    // Writes bufsize() output frames to `output`. Returns the onset value.
    float step(const float *input, float position_pct, float *output);

    void set_stretch_envelope(std::vector<Breakpoint> envelope);
    void clear_stretch_envelope();
    const std::vector<Breakpoint> &stretch_envelope() const;

    void set_onset_detection_sensitivity(float s);

    // Reset internal DSP state (use for seek / loop). Configuration and
    // envelope are preserved.
    void reset();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class OfflineRenderer {
public:
    explicit OfflineRenderer(RenderOptions options = {});

    void set_stretch_envelope(std::vector<Breakpoint> envelope);
    void clear_stretch_envelope();
    const std::vector<Breakpoint> &stretch_envelope() const;

    const RenderOptions &options() const;

    std::vector<float> render_mono(const std::vector<float> &input) const;
    StereoBuffer render_stereo(const std::vector<float> &left,
                               const std::vector<float> &right) const;

    std::size_t estimate_output_frames(std::size_t input_frames) const;

private:
    RenderOptions options_;
    std::vector<Breakpoint> envelope_;
};

} // namespace paulstretch
