#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

#include "paulstretch/paulstretch.h"

// ---------- CLI ----------

static void print_usage(const char *argv0) {
    std::cerr
        << "Usage: " << argv0 << " [options] <input.wav> <output.wav>\n\n"
        << "Options:\n"
        << "  -s, --stretch <ratio>     Stretch ratio (default: 8.0)\n"
        << "  -f, --fft-size <size>     FFT size, must be power of 2 (default: 4096)\n"
        << "  -w, --window <type>       Window: rectangular, hamming, hann,\n"
        << "                            blackman, blackman-harris (default: hann)\n"
        << "  -o, --onset <sensitivity> Onset detection sensitivity 0-1 (default: 0)\n"
        << "  -h, --help                Show this help\n";
}

static paulstretch::Window parse_window(const std::string &s) {
    if (s == "rectangular") return paulstretch::Window::Rectangular;
    if (s == "hamming") return paulstretch::Window::Hamming;
    if (s == "hann") return paulstretch::Window::Hann;
    if (s == "blackman") return paulstretch::Window::Blackman;
    if (s == "blackman-harris") return paulstretch::Window::BlackmanHarris;
    std::cerr << "warning: unknown window '" << s << "', using hann\n";
    return paulstretch::Window::Hann;
}

int main(int argc, char *argv[]) {
    paulstretch::RenderOptions opts;
    std::string input_path, output_path;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-s" || arg == "--stretch") && i + 1 < argc) {
            opts.stretch = std::stof(argv[++i]);
        } else if ((arg == "-f" || arg == "--fft-size") && i + 1 < argc) {
            opts.fft_size = std::stoi(argv[++i]);
        } else if ((arg == "-w" || arg == "--window") && i + 1 < argc) {
            opts.window = parse_window(argv[++i]);
        } else if ((arg == "-o" || arg == "--onset") && i + 1 < argc) {
            opts.onset_detection_sensitivity = std::stof(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return EXIT_SUCCESS;
        } else if (arg[0] == '-') {
            std::cerr << "error: unknown option '" << arg << "'\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        } else if (input_path.empty()) {
            input_path = arg;
        } else if (output_path.empty()) {
            output_path = arg;
        } else {
            std::cerr << "error: unexpected argument '" << arg << "'\n";
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }
    }

    if (input_path.empty() || output_path.empty()) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    // Read input WAV
    drwav wav_in;
    if (!drwav_init_file(&wav_in, input_path.c_str(), nullptr)) {
        std::cerr << "error: cannot open '" << input_path << "'\n";
        return EXIT_FAILURE;
    }

    if (wav_in.channels < 1 || wav_in.channels > 2) {
        std::cerr << "error: only mono/stereo WAV is supported\n";
        drwav_uninit(&wav_in);
        return EXIT_FAILURE;
    }

    auto num_channels = wav_in.channels;
    auto sample_rate = wav_in.sampleRate;
    auto total_frames = wav_in.totalPCMFrameCount;

    std::vector<float> interleaved(total_frames * num_channels);
    drwav_read_pcm_frames_f32(&wav_in, total_frames, interleaved.data());
    drwav_uninit(&wav_in);

    // Deinterleave
    std::vector<float> left(total_frames), right(total_frames);
    if (num_channels == 1) {
        left = std::move(interleaved);
        right = left;
    } else {
        for (drwav_uint64 i = 0; i < total_frames; i++) {
            left[i] = interleaved[i * 2];
            right[i] = interleaved[i * 2 + 1];
        }
    }

    opts.sample_rate = static_cast<float>(sample_rate);

    std::cerr << "Input: " << total_frames << " frames, " << num_channels
              << "ch, " << sample_rate << " Hz\n"
              << "Stretch: " << opts.stretch << "x, FFT: " << opts.fft_size
              << ", Backend: " << paulstretch::fft_backend_name() << "\n";

    paulstretch::OfflineRenderer renderer(opts);

    // Render and write output
    drwav wav_out;
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = num_channels;
    format.sampleRate = sample_rate;
    format.bitsPerSample = 32;

    if (num_channels == 1) {
        auto output = renderer.render_mono(left);
        std::cerr << "Output: " << output.size() << " frames\n";

        if (!drwav_init_file_write(&wav_out, output_path.c_str(), &format, nullptr)) {
            std::cerr << "error: cannot create '" << output_path << "'\n";
            return EXIT_FAILURE;
        }
        drwav_write_pcm_frames(&wav_out, output.size(), output.data());
        drwav_uninit(&wav_out);
    } else {
        auto output = renderer.render_stereo(left, right);
        std::cerr << "Output: " << output.left.size() << " frames\n";

        // Interleave
        auto frames = output.left.size();
        std::vector<float> out_interleaved(frames * 2);
        for (size_t i = 0; i < frames; i++) {
            out_interleaved[i * 2] = output.left[i];
            out_interleaved[i * 2 + 1] = output.right[i];
        }

        if (!drwav_init_file_write(&wav_out, output_path.c_str(), &format, nullptr)) {
            std::cerr << "error: cannot create '" << output_path << "'\n";
            return EXIT_FAILURE;
        }
        drwav_write_pcm_frames(&wav_out, frames, out_interleaved.data());
        drwav_uninit(&wav_out);
    }

    std::cerr << "Wrote " << output_path << "\n";
    return EXIT_SUCCESS;
}
