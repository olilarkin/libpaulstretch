#include <algorithm>
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

	std::size_t estimateOutputFrames(std::size_t input_frames) const {
		return renderer_.estimate_output_frames(input_frames);
	}

private:
	paulstretch::OfflineRenderer renderer_;
};

} // namespace

EMSCRIPTEN_BINDINGS(paulstretch) {
	emscripten::enum_<paulstretch::Window>("Window")
		.value("Rectangular", paulstretch::Window::Rectangular)
		.value("Hamming", paulstretch::Window::Hamming)
		.value("Hann", paulstretch::Window::Hann)
		.value("Blackman", paulstretch::Window::Blackman)
		.value("BlackmanHarris", paulstretch::Window::BlackmanHarris);

	emscripten::class_<WasmOfflineRenderer>("OfflineRenderer")
		.constructor<>()
		.constructor<float, int, float, paulstretch::Window, float>()
		.function("estimateOutputFrames", &WasmOfflineRenderer::estimateOutputFrames)
		.function("renderMono", &WasmOfflineRenderer::renderMono)
		.function("renderStereo", &WasmOfflineRenderer::renderStereo)
		.function("setStretchEnvelope", &WasmOfflineRenderer::setStretchEnvelope)
		.function("clearStretchEnvelope", &WasmOfflineRenderer::clearStretchEnvelope);

	emscripten::function("fftBackendName", &paulstretch::fft_backend_name);
	emscripten::function("fftSimdArch", &paulstretch::fft_simd_arch);
	emscripten::function("fftSimdSize", &paulstretch::fft_simd_size);
}
