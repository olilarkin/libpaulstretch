# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Native build (defaults to PFFFT backend)
cmake -S . -B build
cmake --build build

# Run tests (smoke + render + streaming)
ctest --test-dir build --output-on-failure

# Run a single test
ctest --test-dir build -R paulstretch_streaming_test --output-on-failure

# macOS Accelerate backend
cmake -S . -B build -DPAULSTRETCH_FFT_BACKEND=ACCELERATE

# WebAssembly build (requires Emscripten). Artifacts land directly in npm/dist/
# so `npm pack` in npm/ produces a publishable tarball with no copying.
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
- `PAULSTRETCH_BUILD_TESTS` — build native tests (default: ON)
- `PAULSTRETCH_BUILD_WASM` — build Emscripten target (default: ON; only fires under `emcmake`)
- `PAULSTRETCH_ENABLE_SIMD` — WASM SIMD support, also gates PFFFT NEON on ARM (default: ON)
- `PAULSTRETCH_ENABLE_FAST_MATH` — enable `-ffast-math` (default: OFF)

## Architecture

FFT-based extreme time-stretching library (Paulstretch algorithm). Single translation unit design:

- `include/paulstretch/paulstretch.h` — Public C++20 API. Three primary classes plus options structs and FFT introspection helpers.
- `src/paulstretch.cpp` — Complete implementation (~1.4k lines): FFT backend abstraction (PFFFT/Accelerate/KissFFT), windowing, phase randomization, overlap-add, onset detection, spectral processing (pitch shift, octave mixer, frequency shift, compressor, bandpass/notch filter, harmonics generator, spread, tonal-noise preservation, arbitrary filter), breakpoint envelope interpolation, binaural-beats oscillator, and rendering loop.
- `src/wasm_bindings.cpp` — Emscripten `embind` wrappers (`WasmOfflineRenderer`, `WasmStreamingStretcher`, `WasmBinauralBeatsProcessor`); marshals `Float32Array` and accepts JS option objects in camelCase.
- `examples/cli.cpp` — WAV-to-WAV CLI using dr_wav.

### Public classes

- **`StreamingStretcher`** — Block-based push/pull primitive for realtime use (AudioWorklet / Web Worker). The host calls `next_input_size()` to learn how many input frames `step()` wants next, gathers exactly that many (zero-padding if the source ran out), calls `step()` to produce `bufsize()` output frames, then advances its input cursor by an additional `skip_after_step()` frames. The first `step()` call expects `max_input_chunk()` (= 3 × bufsize) frames for the initial fill; subsequent requests alternate between 0 and `bufsize()` depending on stretch factor and onset detection. `step_without_onset_feedback()` + `apply_onset()` exists for hosts that need to coordinate onsets across channels.
- **`OfflineRenderer`** — Convenience wrapper around `StreamingStretcher` for whole-buffer rendering (`render_mono`, `render_stereo`). Stereo rendering runs two independent `StreamingStretcher`s but synchronizes onset detection (`max` of both channels) so the channels stay aligned.
- **`BinauralBeatsProcessor`** — Independent post-processor that mixes the stretched signal toward mono and adds a sub-audio LFO-style beat between L/R (`set_options`, optional frequency envelope, `process(left, right, nframes, position_pct)`).

### Key data flow

Input samples → window into overlap-add frame → forward FFT → phase randomization → optional spectral processing (`ProcessOptions`: pitch shift, octave mixer, frequency shift, compressor, filter, harmonics, spread, tonal-noise preservation, arbitrary filter) → inverse FFT → overlap-add output buffer → stretched audio. An optional `Breakpoint` envelope (positions 0–1, multiplier values) modulates the stretch ratio per-step based on input cursor position (`position_pct` in `StreamingStretcher::step`).

### Important invariants

- **FFT size rounding:** `fft_size` is auto-rounded to the nearest valid transform size for the active backend (PFFFT has specific size constraints, Accelerate requires power-of-2). After construction, read the sanitized value via `renderer.options().fft_size` / `streaming.options().fft_size`.
- **Hot-set vs. reconfigure:** `set_stretch_factor()` and envelope/process-option setters mutate live DSP state without audible discontinuity. Construction parameters (sample rate, fft_size, window) are immutable — build a new instance to change them. `reset()` clears DSP state for seek/loop but preserves configuration.
- **WASM build output path:** the `paulstretch_wasm` target writes `paulstretch.{js,wasm}` directly into `npm/dist/` (set via `RUNTIME_OUTPUT_DIRECTORY` in `CMakeLists.txt`). Don't rename or move these — `npm/package.json` and the browser smoke test in CI both expect that path.
- **PFFFT** is fetched via CMake FetchContent at a pinned commit (`PAULSTRETCH_PFFFT_GIT_TAG`). KissFFT is bundled in `vendor/kissfft/`. Accelerate links against the system framework.

## CI

`.github/workflows/ci.yml` runs the native test matrix across {ubuntu, macos, windows} × {PFFFT, KISSFFT} (plus ACCELERATE on macOS), and a separate WASM job that builds with SIMD on/off, runs `npm pack`, and executes a headless-Chrome browser smoke test against `tests/browser/index.html`. When modifying workflow files, check status with `gh run list`.

## License

GPLv2
