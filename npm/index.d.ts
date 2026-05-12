// TypeScript types for the Emscripten-generated paulstretch module.
// The .js entry point is an MODULARIZE=1 / EXPORT_ES6=1 factory: call it (with
// optional locateFile) to get a Promise<PaulstretchModule>.

export enum Window {
  Rectangular = 0,
  Hamming = 1,
  Hann = 2,
  Blackman = 3,
  BlackmanHarris = 4,
}

/** Returned by `renderStereo`. */
export interface StereoBuffer {
  left: Float32Array;
  right: Float32Array;
}

export interface OfflineRenderer {
  /** Upper bound on output frame count for the given input length. */
  estimateOutputFrames(inputFrames: number): number;

  renderMono(input: Float32Array): Float32Array;
  renderStereo(left: Float32Array, right: Float32Array): StereoBuffer;

  /**
   * Set a time-varying stretch multiplier. Positions are normalized 0..1
   * over the input duration; values are multipliers on the constructor's
   * `stretch` argument. Breakpoints are sorted internally.
   */
  setStretchEnvelope(positions: Float32Array, values: Float32Array): void;
  clearStretchEnvelope(): void;

  /**
   * Free the underlying C++ object. Failure to call this leaks WASM heap
   * memory — embind objects are not garbage-collected.
   */
  delete(): void;
}

export interface OfflineRendererConstructor {
  new (): OfflineRenderer;
  new (
    stretch: number,
    fftSize: number,
    sampleRate: number,
    window: Window,
    onsetDetectionSensitivity: number,
  ): OfflineRenderer;
}

/**
 * Result of one `StreamingStretcher.step()` call: an output chunk of
 * `bufsize()` frames and the onset detection value.
 */
export interface StreamingStep {
  output: Float32Array;
  onset: number;
}

/**
 * Block-based, push/pull stretching primitive for realtime use (e.g. from a
 * Web Worker driving an AudioWorklet via a ring buffer). The protocol is:
 *
 *   1. Construct with the desired RenderOptions arguments.
 *   2. Call `nextInputSize()` to learn how many input frames `step()` wants.
 *      The first call returns `maxInputChunk()` (initial fill); subsequent
 *      calls return either `0` or `bufsize()`.
 *   3. Gather that many input frames (zero-pad if your source ran out).
 *   4. Call `step(input, positionPct)` where `positionPct` is 0..100 for the
 *      input cursor (used to evaluate the stretch envelope). Receive
 *      `{ output, onset }` where output is a Float32Array of `bufsize()`
 *      frames.
 *   5. Advance your input cursor by `skipAfterStep()` frames after `step()`
 *      (in addition to the frames already consumed).
 */
export interface StreamingStretcher {
  /** Frames the next `step()` call expects: `maxInputChunk()` on first call, then `0` or `bufsize()`. */
  nextInputSize(): number;

  /** Output chunk size. Each `step()` returns this many frames. */
  bufsize(): number;

  /** Largest single input chunk (`3 * bufsize()`); used for the initial fill. */
  maxInputChunk(): number;

  /** Frames to skip in the caller's input cursor after `step()` (in addition to consumed frames). */
  skipAfterStep(): number;

  /**
   * Advance one step.
   * - `input`: Float32Array of `nextInputSize()` frames, or null/undefined if no input is needed.
   * - `positionPct`: input cursor as percent 0..100; used to evaluate the stretch envelope.
   */
  step(input: Float32Array | null | undefined, positionPct: number): StreamingStep;

  /**
   * Advance one step without feeding the returned onset back into this
   * stretcher. Multichannel hosts can call this on every channel, take the
   * maximum onset, then call `applyOnset(maxOnset)` on every channel so their
   * input protocol stays aligned.
   */
  stepWithoutOnsetFeedback(input: Float32Array | null | undefined, positionPct: number): StreamingStep;
  applyOnset(onset: number): void;

  setStretchEnvelope(positions: Float32Array, values: Float32Array): void;
  clearStretchEnvelope(): void;

  /** Hot-swap the base stretch factor without resetting DSP state. */
  setStretchFactor(stretch: number): void;

  setOnsetDetectionSensitivity(s: number): void;

  /** Reset internal DSP state (seek/loop). Configuration and envelope are preserved. */
  reset(): void;

  /** Free the underlying C++ object. Failure to call this leaks WASM heap memory. */
  delete(): void;
}

export interface StreamingStretcherConstructor {
  new (
    stretch: number,
    fftSize: number,
    sampleRate: number,
    window: Window,
    onsetDetectionSensitivity: number,
  ): StreamingStretcher;
}

export interface PaulstretchModule {
  OfflineRenderer: OfflineRendererConstructor;
  StreamingStretcher: StreamingStretcherConstructor;
  Window: typeof Window;
  fftBackendName(): string;
  fftSimdArch(): string;
  fftSimdSize(): number;
}

export interface ModuleFactoryOptions {
  /**
   * Override where the runtime looks for paulstretch.wasm. Useful when the
   * .wasm is served from a different path than the .js, or when bundling
   * for the browser.
   */
  locateFile?: (path: string, scriptDirectory: string) => string;

  /** Emscripten-standard hooks. */
  print?: (text: string) => void;
  printErr?: (text: string) => void;
}

declare const PaulstretchModuleFactory: (
  options?: ModuleFactoryOptions,
) => Promise<PaulstretchModule>;

export default PaulstretchModuleFactory;
