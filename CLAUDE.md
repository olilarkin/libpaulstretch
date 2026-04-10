# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Native build (defaults to PFFFT backend)
cmake -S . -B build
cmake --build build

# Run tests
ctest --test-dir build --output-on-failure

# macOS Accelerate backend
cmake -S . -B build -DPAULSTRETCH_FFT_BACKEND=ACCELERATE

# WebAssembly build (requires Emscripten)
emcmake cmake -S . -B build-wasm
cmake --build build-wasm

# WASM without SIMD
emcmake cmake -S . -B build-wasm -DPAULSTRETCH_ENABLE_SIMD=OFF
```

The CLI example reads/writes WAV files via dr_wav (header-only, in `vendor/dr_libs/`):
```bash
./build/paulstretch -s 8 -f 4096 -w hann input.wav output.wav
```

### CMake Options

- `PAULSTRETCH_FFT_BACKEND` — `PFFFT` (default, SIMD-capable), `ACCELERATE` (macOS/iOS vDSP), or `KISSFFT` (bundled fallback)
- `PAULSTRETCH_BUILD_EXAMPLES` — build CLI example (default: ON)
- `PAULSTRETCH_BUILD_TESTS` — build smoke test (default: ON)
- `PAULSTRETCH_BUILD_WASM` — build Emscripten target (default: ON)
- `PAULSTRETCH_ENABLE_SIMD` — WASM SIMD support (default: ON)
- `PAULSTRETCH_ENABLE_FAST_MATH` — enable `-ffast-math` (default: OFF)

## Architecture

FFT-based extreme time-stretching library (Paulstretch algorithm). Single translation unit design:

- `include/paulstretch/paulstretch.h` — Public C++20 API: `OfflineRenderer`, `RenderOptions`, `Window`, `Breakpoint`, `StereoBuffer`.
- `src/paulstretch.cpp` — Complete implementation: FFT backend abstraction (PFFFT/Accelerate/KissFFT), core stretch algorithm (windowing, phase randomization, overlap-add, onset detection), breakpoint envelope interpolation, and rendering loop.
- `src/wasm_bindings.cpp` — Emscripten bindings exposing `OfflineRenderer` to JavaScript via `embind`.
- `examples/cli.cpp` — WAV-to-WAV CLI using dr_wav.

**Key data flow:** Input samples → FFT analysis with windowing → phase randomization → inverse FFT → overlap-add output buffer → stretched audio. An optional breakpoint envelope modulates the stretch ratio per-chunk based on input position.

**FFT size rounding:** `fft_size` is auto-rounded to the nearest valid transform size for the active backend (PFFFT has specific size constraints, Accelerate requires power-of-2). The sanitized value is accessible via `renderer.options().fft_size`.

**Stereo rendering** processes left/right channels through independent `Stretcher` instances but synchronizes onset detection (max of both channels) to keep channels aligned.

**PFFFT** is fetched via CMake FetchContent (pinned commit). KissFFT is bundled in `vendor/kissfft/`.

## License

GPLv2
