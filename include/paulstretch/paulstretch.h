#pragma once

#include <cstddef>
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
