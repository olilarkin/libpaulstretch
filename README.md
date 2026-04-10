# libpaulstretch

A portable C++20 library implementing the [Paulstretch](http://hypermammut.sourceforge.net/paulstretch/) extreme time-stretching algorithm. Includes an Emscripten/WebAssembly target.

## Native build

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

### FFT backends

The library supports three FFT backends, selected via the `PAULSTRETCH_FFT_BACKEND` CMake option:

| Backend | Value | Description |
|---|---|---|
| **PFFFT** | `PFFFT` (default) | Fetched via CMake FetchContent from [marton78/pffft](https://github.com/marton78/pffft). SIMD-capable (SSE, NEON, WASM SIMD). Has specific transform-size constraints — the renderer auto-rounds `fft_size` to the nearest valid size. |
| **Accelerate** | `ACCELERATE` | Uses Apple's vDSP via the Accelerate framework. macOS/iOS only. Requires power-of-2 FFT sizes (auto-rounded). |
| **KissFFT** | `KISSFFT` | Bundled in `vendor/kissfft/`. Pure C, no SIMD — useful as a portable fallback when neither PFFFT nor Accelerate is available. |

```bash
# Default (PFFFT)
cmake -S . -B build

# macOS Accelerate
cmake -S . -B build -DPAULSTRETCH_FFT_BACKEND=ACCELERATE

# KissFFT fallback
cmake -S . -B build -DPAULSTRETCH_FFT_BACKEND=KISSFFT
```

## WebAssembly build

```bash
emcmake cmake -S . -B build-wasm
cmake --build build-wasm
```

Output files: `build-wasm/dist/paulstretch.js` and `build-wasm/dist/paulstretch.wasm`.

SIMD is enabled by default (`-msimd128`). Disable with:

```bash
emcmake cmake -S . -B build-wasm -DPAULSTRETCH_ENABLE_SIMD=OFF
```

## C++ usage

```cpp
#include "paulstretch/paulstretch.h"

paulstretch::OfflineRenderer renderer({
    .stretch = 8.0f,
    .fft_size = 4096,
    .sample_rate = 48000.0f,
    .window = paulstretch::Window::Hann,
});

std::vector<float> output = renderer.render_mono(input);
```

### Stretch envelope

Apply a time-varying stretch multiplier with breakpoints. Each breakpoint has a `position` (0–1 normalized time) and a `value` (stretch multiplier relative to the base stretch factor):

```cpp
renderer.set_stretch_envelope({
    {0.0f, 1.0f},   // start: 1x base stretch
    {0.5f, 4.0f},   // middle: 4x base stretch
    {1.0f, 1.0f},   // end: 1x base stretch
});

std::vector<float> output = renderer.render_mono(input);
renderer.clear_stretch_envelope(); // revert to uniform
```

### FFT backend introspection

```cpp
std::string backend = paulstretch::fft_backend_name(); // "PFFFT"
std::string simd = paulstretch::fft_simd_arch();       // "NEON", "SSE", ...
int width = paulstretch::fft_simd_size();              // 4
```

## Node.js / WASM usage

```js
import createPaulstretchModule from "./build-wasm/dist/paulstretch.js";

const Module = await createPaulstretchModule();
const renderer = new Module.OfflineRenderer(8.0, 4096, 48000, Module.Window.Hann, 0.0);

const output = renderer.renderMono(input);
const stereo = renderer.renderStereo(leftChannel, rightChannel);
```

### Stretch envelope

Pass parallel arrays of positions (0–1) and multiplier values:

```js
renderer.setStretchEnvelope([0, 0.5, 1.0], [1.0, 4.0, 1.0]);
const output = renderer.renderMono(input);
renderer.clearStretchEnvelope();
```

## Notes

- PFFFT only supports specific transform sizes. The renderer rounds `fft_size` to the nearest valid size automatically.
- GPLv2-licensed.
