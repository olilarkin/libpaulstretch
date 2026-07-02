# libpaulstretch

A portable C++20 library implementing the [Paulstretch](http://hypermammut.sourceforge.net/paulstretch/) extreme time-stretching algorithm. Ships an offline renderer, a realtime streaming primitive, optional spectral processing (pitch shift, octave mixer, frequency shift, compressor, filter, harmonics, spread, tonal-noise preservation, arbitrary filter), a binaural-beats post-processor, and an Emscripten/WebAssembly target.

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

Output files land directly in `npm/dist/paulstretch.js` and `npm/dist/paulstretch.wasm` so that `cd npm && npm pack` produces a publishable tarball with no extra copying.

SIMD is enabled by default (`-msimd128`). Disable with:

```bash
emcmake cmake -S . -B build-wasm -DPAULSTRETCH_ENABLE_SIMD=OFF
```

## C++ usage

### Offline rendering

```cpp
#include "paulstretch/paulstretch.h"

paulstretch::OfflineRenderer renderer({
    .stretch = 8.0f,
    .fft_size = 4096,
    .sample_rate = 48000.0f,
    .window = paulstretch::Window::Hann,
});

std::vector<float> output = renderer.render_mono(input);
auto [left, right] = renderer.render_stereo(left_in, right_in);
```

Stereo rendering runs two independent stretchers but synchronizes onset detection across channels so they stay phase-aligned.

#### Very long outputs (chunked rendering)

`render_mono` returns the whole result in one `std::vector`. For an extreme stretch — a few seconds blown up several hundred times into an hour-plus of audio — that buffer can be enormous, and on the WebAssembly build it has to live in linear memory twice over (the C++ vector plus the returned `Float32Array`), which can exceed the WASM heap and abort. `render_mono_chunked` / `render_stereo_chunked` run the identical algorithm but hand each `bufsize()`-frame chunk to a sink as it is produced, so peak memory stays bounded regardless of output length:

```cpp
renderer.render_mono_chunked(input, [&](const float *data, int frames) {
    // Consume the chunk — append to a buffer, write to disk, feed an encoder.
});

renderer.render_stereo_chunked(left_in, right_in,
    [&](const float *left, const float *right, int frames) { /* ... */ });
```

### Streaming (realtime) rendering

`StreamingStretcher` is a block-based push/pull primitive for realtime hosts (audio callback, AudioWorklet, Web Worker). The host gathers exactly the number of input frames the stretcher asks for, calls `step()` to produce one output chunk, then advances its input cursor by the additional skip distance:

```cpp
paulstretch::StreamingStretcher s({
    .stretch = 8.0f,
    .fft_size = 4096,
    .sample_rate = 48000.0f,
});

std::vector<float> out(s.bufsize());

while (rendering) {
    const int want = s.next_input_size();    // first call: 3 * bufsize() for initial fill
    std::vector<float> in(want);
    read_from_source(in.data(), want);       // zero-pad if source ran out

    const float position_pct = 100.0f * input_cursor / total_input_frames;
    s.step(in.data(), position_pct, out.data());
    write_to_output(out.data(), s.bufsize());

    input_cursor += want + s.skip_after_step();
}
```

`set_stretch_factor()` hot-swaps the base stretch ratio without resetting DSP state, and `reset()` clears state for seek/loop while preserving configuration.

### Stretch envelope

Apply a time-varying stretch multiplier with breakpoints. Each breakpoint has a `position` (0–1 normalized time) and a `value` (multiplier on the base stretch factor):

```cpp
renderer.set_stretch_envelope({
    {0.0f, 1.0f},   // start: 1× base stretch
    {0.5f, 4.0f},   // middle: 4× base stretch
    {1.0f, 1.0f},   // end: 1× base stretch
});

std::vector<float> output = renderer.render_mono(input);
renderer.clear_stretch_envelope(); // revert to uniform
```

### Spectral processing

`ProcessOptions` enables effects that run on the FFT bins between phase randomization and the inverse transform. All effects default to disabled; set the matching `*_enabled` flag plus its parameters:

```cpp
paulstretch::ProcessOptions p;
p.pitch_shift_enabled = true;
p.pitch_shift_cents = 700;             // up a perfect fifth

p.filter_enabled = true;
p.filter_low_hz = 200.0f;
p.filter_high_hz = 4000.0f;

p.harmonics_enabled = true;
p.harmonics_frequency_hz = 110.0f;
p.harmonics_count = 8;

renderer.set_process_options(p);
```

Available effects: pitch shift, octave mixer (–2/–1/0/+1/+1.5/+2), frequency shift, compressor, bandpass/notch filter, harmonics generator, stereo spread, tonal-noise preservation, and an arbitrary breakpoint-shaped filter (`set_arbitrary_filter` + `arbitrary_filter_enabled`). See `include/paulstretch/paulstretch.h` for the full struct.

### Binaural beats

Post-process the stretched output to add a sub-audio beat between L/R channels:

```cpp
paulstretch::BinauralBeatsProcessor bb(48000.0f);
bb.set_options({
    .enabled = true,
    .stereo_mode = paulstretch::BinauralStereoMode::LeftRight,
    .mono = 0.5f,              // mix toward mono before applying the beat
    .beat_frequency_hz = 8.0f, // alpha range
});

bb.process(left.data(), right.data(), nframes, position_pct);
```

The beat frequency can itself be automated with `set_frequency_envelope()`.

### FFT backend introspection

```cpp
std::string backend = paulstretch::fft_backend_name(); // "PFFFT"
std::string simd = paulstretch::fft_simd_arch();       // "NEON", "SSE", ...
int width = paulstretch::fft_simd_size();              // 4
```

## Node.js / WASM usage

The Emscripten build is published as [`@olilarkin/paulstretch-wasm`](https://github.com/olilarkin/libpaulstretch/pkgs/npm/paulstretch-wasm). See [`npm/README.md`](npm/README.md) for installation, the full JS API (offline, streaming, envelope, spectral processing, binaural beats), and bundler notes. Type definitions live in [`npm/index.d.ts`](npm/index.d.ts).

The C++ and JS APIs are 1:1 — JS methods use camelCase versions of the C++ names, and `setProcessOptions` accepts a plain object instead of a `ProcessOptions` struct.

## Notes

- PFFFT only supports specific transform sizes. The renderer rounds `fft_size` to the nearest valid size automatically; read the effective value from `renderer.options().fft_size` (C++) after construction.
- Sample rate, FFT size, and window are immutable after construction — build a new renderer to change them. Stretch factor, envelopes, and process options are hot-settable.
- GPLv2-licensed.
