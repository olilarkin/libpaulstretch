#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "paulstretch/paulstretch.h"

namespace {

std::vector<float> make_sine(std::size_t frames, float sample_rate, float frequency) {
	std::vector<float> result(frames, 0.0f);
	const float tau = 6.28318530717958647692f;
	for (std::size_t i = 0; i < frames; i++) {
		float phase = tau * frequency * static_cast<float>(i) / sample_rate;
		result[i] = 0.5f * std::sin(phase);
	}
	return result;
}

bool all_finite(const std::vector<float> &samples) {
	for (float s : samples) {
		if (!std::isfinite(s)) return false;
	}
	return true;
}

} // namespace

int main() {
	std::cout << "backend=" << paulstretch::fft_backend_name()
	          << " simd=" << paulstretch::fft_simd_arch()
	          << " width=" << paulstretch::fft_simd_size() << '\n';

	const paulstretch::RenderOptions options{
		.stretch = 6.0f,
		.fft_size = 1024,
		.sample_rate = 48000.0f,
		.window = paulstretch::Window::Hann,
		.onset_detection_sensitivity = 0.0f,
	};

	paulstretch::OfflineRenderer renderer(options);
	const auto mono_input = make_sine(4800, options.sample_rate, 220.0f);
	const auto mono_output = renderer.render_mono(mono_input);

	if (mono_output.empty()) {
		std::cerr << "render_mono returned no samples\n";
		return EXIT_FAILURE;
	}
	if (mono_output.size() <= mono_input.size()) {
		std::cerr << "render_mono did not stretch the signal\n";
		return EXIT_FAILURE;
	}
	if (!all_finite(mono_output)) {
		std::cerr << "render_mono produced non-finite samples\n";
		return EXIT_FAILURE;
	}

	const auto stereo_output = renderer.render_stereo(mono_input, mono_input);
	if (stereo_output.left.size() != stereo_output.right.size()) {
		std::cerr << "stereo channel lengths do not match\n";
		return EXIT_FAILURE;
	}
	if (stereo_output.left.empty()) {
		std::cerr << "render_stereo returned no samples\n";
		return EXIT_FAILURE;
	}
	if (!all_finite(stereo_output.left) || !all_finite(stereo_output.right)) {
		std::cerr << "render_stereo produced non-finite samples\n";
		return EXIT_FAILURE;
	}

	std::cout << "mono_frames=" << mono_output.size()
	          << " stereo_frames=" << stereo_output.left.size() << '\n';
	return EXIT_SUCCESS;
}
