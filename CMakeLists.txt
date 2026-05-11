cmake_minimum_required(VERSION 3.20)

project(paulstretch_cpp LANGUAGES C CXX)

include(FetchContent)

option(PAULSTRETCH_BUILD_EXAMPLES "Build CLI example" ON)
option(PAULSTRETCH_BUILD_TESTS "Build the native smoke test" ON)
option(PAULSTRETCH_BUILD_WASM "Build the Emscripten Node.js/WebAssembly module" ON)
option(PAULSTRETCH_ENABLE_SIMD "Enable SIMD for the Emscripten build" ON)
option(PAULSTRETCH_ENABLE_FAST_MATH "Enable fast-math for the modern targets" OFF)
set(PAULSTRETCH_FFT_BACKEND "PFFFT" CACHE STRING "FFT backend for the portable core")
set_property(CACHE PAULSTRETCH_FFT_BACKEND PROPERTY STRINGS PFFFT KISSFFT ACCELERATE)
set(PAULSTRETCH_PFFFT_GIT_TAG "0b4ee12c4ba45a4a8e567550c16d96d1679f50ce" CACHE STRING "Pinned marton78/pffft revision")

add_library(paulstretch_core STATIC
	src/paulstretch.cpp
)

target_compile_features(paulstretch_core PUBLIC cxx_std_20)
target_include_directories(paulstretch_core
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
)
set_target_properties(paulstretch_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

if (PAULSTRETCH_FFT_BACKEND STREQUAL "PFFFT")
	set(PFFFT_USE_TYPE_FLOAT ON CACHE BOOL "" FORCE)
	set(PFFFT_USE_TYPE_DOUBLE OFF CACHE BOOL "" FORCE)
	set(PFFFT_USE_SIMD ${PAULSTRETCH_ENABLE_SIMD} CACHE BOOL "" FORCE)
	set(PFFFT_USE_SCALAR_VECT ON CACHE BOOL "" FORCE)
	set(PFFFT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
	set(PFFFT_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
	set(PFFFT_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
	set(INSTALL_PFFFT OFF CACHE BOOL "" FORCE)
	set(INSTALL_PFDSP OFF CACHE BOOL "" FORCE)
	set(INSTALL_PFFASTCONV OFF CACHE BOOL "" FORCE)

	FetchContent_Declare(pffft
		GIT_REPOSITORY https://github.com/marton78/pffft.git
		GIT_TAG ${PAULSTRETCH_PFFFT_GIT_TAG}
	)
	FetchContent_MakeAvailable(pffft)

	if (PAULSTRETCH_ENABLE_SIMD AND (EMSCRIPTEN OR CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|arm64|aarch64)"))
		target_compile_definitions(PFFFT PRIVATE PFFFT_ENABLE_NEON=1)
	endif()
	if (EMSCRIPTEN AND PAULSTRETCH_ENABLE_SIMD)
		target_compile_options(PFFFT PRIVATE -msimd128)
	endif()

	target_link_libraries(paulstretch_core PUBLIC PFFFT::PFFFT)
	target_compile_definitions(paulstretch_core PUBLIC PAULSTRETCH_USE_PFFFT=1)
elseif (PAULSTRETCH_FFT_BACKEND STREQUAL "ACCELERATE")
	if (NOT APPLE)
		message(FATAL_ERROR "ACCELERATE backend requires macOS/iOS (Apple Accelerate framework)")
	endif()
	find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
	target_link_libraries(paulstretch_core PUBLIC ${ACCELERATE_FRAMEWORK})
	target_compile_definitions(paulstretch_core PUBLIC PAULSTRETCH_USE_ACCELERATE=1)
elseif (PAULSTRETCH_FFT_BACKEND STREQUAL "KISSFFT")
	target_sources(paulstretch_core PRIVATE
		vendor/kissfft/kiss_fft.c
		vendor/kissfft/kiss_fftr.c
	)
	target_compile_definitions(paulstretch_core PUBLIC KISSFFT)
	target_include_directories(paulstretch_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/vendor/kissfft)
else()
	message(FATAL_ERROR "Unsupported PAULSTRETCH_FFT_BACKEND='${PAULSTRETCH_FFT_BACKEND}'")
endif()

if (PAULSTRETCH_ENABLE_FAST_MATH)
	target_compile_options(paulstretch_core PRIVATE
		$<$<COMPILE_LANGUAGE:CXX>:-ffast-math>
		$<$<COMPILE_LANGUAGE:C>:-ffast-math>
	)
endif()

if (PAULSTRETCH_BUILD_TESTS AND NOT EMSCRIPTEN)
	enable_testing()
	add_executable(paulstretch_smoke tests/smoke_test.cpp)
	target_link_libraries(paulstretch_smoke PRIVATE paulstretch_core)
	add_test(NAME paulstretch_smoke COMMAND paulstretch_smoke)

	add_executable(paulstretch_render_test tests/render_test.cpp)
	target_link_libraries(paulstretch_render_test PRIVATE paulstretch_core)
	target_include_directories(paulstretch_render_test PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/vendor/dr_libs
	)
	add_test(NAME paulstretch_render_test
		COMMAND paulstretch_render_test
			${CMAKE_CURRENT_SOURCE_DIR}/tests/plenty_of_unknown.wav
	)

	add_executable(paulstretch_streaming_test tests/streaming_test.cpp)
	target_link_libraries(paulstretch_streaming_test PRIVATE paulstretch_core)
	target_include_directories(paulstretch_streaming_test PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/vendor/dr_libs
	)
	add_test(NAME paulstretch_streaming_test
		COMMAND paulstretch_streaming_test
			${CMAKE_CURRENT_SOURCE_DIR}/tests/plenty_of_unknown.wav
	)
endif()

if (PAULSTRETCH_BUILD_EXAMPLES AND NOT EMSCRIPTEN)
	add_executable(paulstretch examples/cli.cpp)
	target_link_libraries(paulstretch PRIVATE paulstretch_core)
	target_include_directories(paulstretch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/vendor/dr_libs)
endif()

if (EMSCRIPTEN AND PAULSTRETCH_BUILD_WASM)
	add_executable(paulstretch_wasm src/wasm_bindings.cpp)
	target_link_libraries(paulstretch_wasm PRIVATE paulstretch_core)
	target_compile_features(paulstretch_wasm PRIVATE cxx_std_20)
	set_target_properties(paulstretch_wasm PROPERTIES
		OUTPUT_NAME "paulstretch"
		SUFFIX ".js"
		# Build straight into the npm package's dist/ so `npm pack` after a
		# build produces a publishable tarball with no extra copying.
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/npm/dist"
	)

	target_link_options(paulstretch_wasm PRIVATE
		-O3
		--bind
		-sALLOW_MEMORY_GROWTH=1
		# Default ENVIRONMENT (web,webview,worker,node) — works in browsers,
		# bundlers (Vite/Webpack), Node, and Web Workers.
		-sEXPORT_ES6=1
		-sFILESYSTEM=0
		-sMODULARIZE=1
		-sWASM=1
	)

	if (PAULSTRETCH_ENABLE_SIMD)
		target_compile_options(paulstretch_core PRIVATE -msimd128)
		target_compile_options(paulstretch_wasm PRIVATE -msimd128)
		target_link_options(paulstretch_wasm PRIVATE -msimd128)
	endif()

	if (PAULSTRETCH_ENABLE_FAST_MATH)
		target_link_options(paulstretch_wasm PRIVATE -ffast-math)
	endif()
endif()
