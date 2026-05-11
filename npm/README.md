# paulstretch-wasm

[Paulstretch](https://github.com/paulnasca/paulstretch_cpp) extreme time-stretching, compiled to WebAssembly. Works in Node, modern browsers, and Web Workers.

## Install

```bash
npm install paulstretch-wasm
```

## Usage

```js
import PaulstretchModule from 'paulstretch-wasm';

const Module = await PaulstretchModule();

const stretch = 8;
const fftSize = 4096;
const sampleRate = 44100;
const renderer = new Module.OfflineRenderer(
  stretch,
  fftSize,
  sampleRate,
  Module.Window.Hann,
  0, // onset detection sensitivity, 0..1
);

const input = new Float32Array(/* your audio, [-1, 1] floats */);
const output = renderer.renderMono(input);

// Always free embind objects — they don't get GC'd automatically.
renderer.delete();
```

### Stereo

```js
const { left, right } = renderer.renderStereo(leftIn, rightIn);
```

### Time-varying stretch (breakpoint envelope)

Positions are normalized `0..1` over the input. Values multiply the `stretch` you passed to the constructor.

```js
renderer.setStretchEnvelope(
  new Float32Array([0.0, 0.5, 1.0]),
  new Float32Array([1.0, 4.0, 1.0]), // 1× → 4× → 1× over the input
);
const out = renderer.renderMono(input);
```

### Browser bundlers

When bundling with Vite/Webpack/esbuild, the runtime may need help finding `paulstretch.wasm`:

```js
import wasmUrl from 'paulstretch-wasm/paulstretch.wasm?url';
import PaulstretchModule from 'paulstretch-wasm';

const Module = await PaulstretchModule({
  locateFile: (path) => (path.endsWith('.wasm') ? wasmUrl : path),
});
```

## Building from source

```bash
emcmake cmake -S . -B build-wasm
cmake --build build-wasm
# outputs land in npm/dist/
cd npm && npm pack
```

## License

GPL-2.0 — same as the upstream Paulstretch algorithm.
