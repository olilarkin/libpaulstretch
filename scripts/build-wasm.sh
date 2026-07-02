#!/usr/bin/env bash
#
# Build BOTH WebAssembly artifacts and assemble the npm package's dist/:
#
#   npm/dist/paulstretch.js           - Emscripten glue (identical for both wasm)
#   npm/dist/paulstretch.wasm         - SIMD build (-msimd128, WASM_SIMD128 PFFFT)
#   npm/dist/paulstretch.nosimd.wasm  - scalar build (no SIMD)
#
# The scalar wasm is a fallback for WebViews without WASM SIMD support (e.g.
# macOS WKWebView before Safari 16.4 / macOS 13). Consumers feature-detect SIMD
# at runtime and load the matching .wasm via Emscripten's `locateFile`. Because
# both wasm are built with the SAME Emscripten version, the generated glue is
# byte-identical, so one glue drives either binary.
#
# Requires emcc/emcmake on PATH (run `source <emsdk>/emsdk_env.sh` first, or use
# the setup-emsdk CI action).
#
# Usage: scripts/build-wasm.sh

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

COMMON_ARGS=(
  -DCMAKE_BUILD_TYPE=Release
  -DPAULSTRETCH_BUILD_WASM=ON
  -DPAULSTRETCH_BUILD_TESTS=OFF
  -DPAULSTRETCH_BUILD_EXAMPLES=OFF
)

# Both variants link straight into npm/dist as paulstretch.{js,wasm}. We delete
# that output before each build so an up-to-date build tree still refreshes it —
# otherwise CMake skips the relink and we'd copy a stale wasm from the other
# variant. Deleting only the final artifact forces a cheap relink, not a full
# recompile.

# 1. Scalar (non-SIMD) build first, then stash it before the SIMD build
#    overwrites paulstretch.wasm.
echo "==> Configuring scalar (SIMD=OFF) build"
emcmake cmake -S . -B build-wasm-nosimd "${COMMON_ARGS[@]}" -DPAULSTRETCH_ENABLE_SIMD=OFF
echo "==> Building scalar wasm"
rm -f npm/dist/paulstretch.js npm/dist/paulstretch.wasm
cmake --build build-wasm-nosimd --config Release --parallel
cp npm/dist/paulstretch.wasm npm/dist/paulstretch.nosimd.wasm

# 2. SIMD build. This leaves npm/dist/paulstretch.{js,wasm} as the SIMD pair,
#    which is the canonical glue + primary artifact.
echo "==> Configuring SIMD (SIMD=ON) build"
emcmake cmake -S . -B build-wasm "${COMMON_ARGS[@]}" -DPAULSTRETCH_ENABLE_SIMD=ON
echo "==> Building SIMD wasm"
rm -f npm/dist/paulstretch.js npm/dist/paulstretch.wasm
cmake --build build-wasm --config Release --parallel

# Sanity: the two binaries must actually differ (SIMD vs scalar). Identical
# bytes means the stale-output hazard bit us again.
if cmp -s npm/dist/paulstretch.wasm npm/dist/paulstretch.nosimd.wasm; then
  echo "error: paulstretch.wasm and paulstretch.nosimd.wasm are identical" >&2
  echo "       (scalar build did not produce a distinct binary)" >&2
  exit 1
fi

echo "==> Assembled npm/dist:"
ls -la npm/dist/paulstretch.js npm/dist/paulstretch.wasm npm/dist/paulstretch.nosimd.wasm
