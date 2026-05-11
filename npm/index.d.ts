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

export interface PaulstretchModule {
  OfflineRenderer: OfflineRendererConstructor;
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
