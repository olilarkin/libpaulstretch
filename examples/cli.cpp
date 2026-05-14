#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

#include "paulstretch/paulstretch.h"

using paulstretch::Breakpoint;

// ---------- parsing helpers ----------

static std::vector<float> split_floats(const std::string &s) {
    std::vector<float> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::stof(tok));
    }
    return out;
}

static std::vector<Breakpoint> parse_breakpoints(const std::string &s) {
    std::vector<Breakpoint> bps;
    std::stringstream ss(s);
    std::string pair;
    while (std::getline(ss, pair, ';')) {
        if (pair.empty()) continue;
        auto comma = pair.find(',');
        if (comma == std::string::npos) {
            throw std::runtime_error("breakpoint missing comma: '" + pair + "'");
        }
        Breakpoint bp;
        bp.position = std::stof(pair.substr(0, comma));
        bp.value = std::stof(pair.substr(comma + 1));
        bps.push_back(bp);
    }
    return bps;
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

static paulstretch::BinauralStereoMode parse_binaural_mode(const std::string &s) {
    if (s == "lr") return paulstretch::BinauralStereoMode::LeftRight;
    if (s == "rl") return paulstretch::BinauralStereoMode::RightLeft;
    if (s == "sym") return paulstretch::BinauralStereoMode::Symmetric;
    std::cerr << "warning: unknown binaural mode '" << s << "', using lr\n";
    return paulstretch::BinauralStereoMode::LeftRight;
}

static void print_usage(const char *argv0) {
    std::cerr <<
        "Usage: " << argv0 << " [options] <input.wav> <output.wav>\n\n"
        "Basic:\n"
        "  -s, --stretch <ratio>           Stretch ratio (default: 8.0)\n"
        "  -f, --fft-size <size>           FFT size (default: 4096)\n"
        "  -w, --window <type>             rectangular|hamming|hann|blackman|blackman-harris\n"
        "  -o, --onset <0-1>               Onset detection sensitivity (default: 0)\n"
        "      --stretch-envelope <bp>     Per-position stretch multipliers,\n"
        "                                  format \"pos,val;pos,val;...\" (pos 0..1)\n\n"
        "Spectral processing:\n"
        "      --pitch-shift <cents>       Pitch shift in cents\n"
        "      --octave <m2,m1,0,p1,p1.5,p2>\n"
        "                                  Octave mixer levels (6 floats)\n"
        "      --frequency-shift <hz>      Linear frequency shift\n"
        "      --compressor <power>        Spectral compressor power\n"
        "      --filter <low,high[,damp]>  Bandpass filter (Hz, optional high damp)\n"
        "      --filter-stop               Make --filter a notch instead of bandpass\n"
        "      --harmonics <freq,bw_cents,count[,gauss]>\n"
        "      --spread <bandwidth>        Spectral spread bandwidth\n"
        "      --tonal-noise <preserve,bw> Tonal/noise preservation\n"
        "      --arbitrary-filter <bp>     Frequency curve breakpoints\n"
        "                                  (pos = normalized bin 0..1, val = gain)\n\n"
        "Binaural beats (post-process; output is forced to stereo):\n"
        "      --binaural <hz>             Beat frequency in Hz (enables binaural)\n"
        "      --binaural-mono <0-1>       Mono mix amount (default 0.5)\n"
        "      --binaural-mode <lr|rl|sym>\n"
        "      --binaural-envelope <bp>    Beat frequency envelope (pos 0..1, val = Hz)\n\n"
        "  -h, --help                      Show this help\n";
}

// ---------- main ----------

int main(int argc, char *argv[]) {
    paulstretch::RenderOptions opts;
    paulstretch::ProcessOptions proc;
    paulstretch::BinauralBeatsOptions binaural;
    std::vector<Breakpoint> stretch_env;
    std::vector<Breakpoint> arb_filter;
    std::vector<Breakpoint> binaural_env;
    std::string input_path, output_path;

    auto need_value = [&](int &i, const std::string &flag) -> std::string {
        if (i + 1 >= argc) {
            throw std::runtime_error(flag + " requires a value");
        }
        return std::string(argv[++i]);
    };

    try {
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "-s" || arg == "--stretch") {
                opts.stretch = std::stof(need_value(i, arg));
            } else if (arg == "-f" || arg == "--fft-size") {
                opts.fft_size = std::stoi(need_value(i, arg));
            } else if (arg == "-w" || arg == "--window") {
                opts.window = parse_window(need_value(i, arg));
            } else if (arg == "-o" || arg == "--onset") {
                opts.onset_detection_sensitivity = std::stof(need_value(i, arg));
            } else if (arg == "--stretch-envelope") {
                stretch_env = parse_breakpoints(need_value(i, arg));
            } else if (arg == "--pitch-shift") {
                proc.pitch_shift_enabled = true;
                proc.pitch_shift_cents = std::stoi(need_value(i, arg));
            } else if (arg == "--octave") {
                auto v = split_floats(need_value(i, arg));
                if (v.size() != 6) {
                    throw std::runtime_error("--octave needs 6 comma-separated values");
                }
                proc.octave_enabled = true;
                proc.octave_minus2 = v[0];
                proc.octave_minus1 = v[1];
                proc.octave_0 = v[2];
                proc.octave_plus1 = v[3];
                proc.octave_plus15 = v[4];
                proc.octave_plus2 = v[5];
            } else if (arg == "--frequency-shift") {
                proc.frequency_shift_enabled = true;
                proc.frequency_shift_hz = std::stoi(need_value(i, arg));
            } else if (arg == "--compressor") {
                proc.compressor_enabled = true;
                proc.compressor_power = std::stof(need_value(i, arg));
            } else if (arg == "--filter") {
                auto v = split_floats(need_value(i, arg));
                if (v.size() < 2 || v.size() > 3) {
                    throw std::runtime_error("--filter needs low,high[,damp]");
                }
                proc.filter_enabled = true;
                proc.filter_low_hz = v[0];
                proc.filter_high_hz = v[1];
                if (v.size() == 3) proc.filter_high_damp = v[2];
            } else if (arg == "--filter-stop") {
                proc.filter_stop = true;
            } else if (arg == "--harmonics") {
                auto v = split_floats(need_value(i, arg));
                if (v.size() < 3 || v.size() > 4) {
                    throw std::runtime_error("--harmonics needs freq,bw_cents,count[,gauss]");
                }
                proc.harmonics_enabled = true;
                proc.harmonics_frequency_hz = v[0];
                proc.harmonics_bandwidth_cents = v[1];
                proc.harmonics_count = static_cast<int>(v[2]);
                if (v.size() == 4) proc.harmonics_gauss = v[3] != 0.0f;
            } else if (arg == "--spread") {
                proc.spread_enabled = true;
                proc.spread_bandwidth = std::stof(need_value(i, arg));
            } else if (arg == "--tonal-noise") {
                auto v = split_floats(need_value(i, arg));
                if (v.size() != 2) {
                    throw std::runtime_error("--tonal-noise needs preserve,bandwidth");
                }
                proc.tonal_noise_enabled = true;
                proc.tonal_noise_preserve = v[0];
                proc.tonal_noise_bandwidth = v[1];
            } else if (arg == "--arbitrary-filter") {
                arb_filter = parse_breakpoints(need_value(i, arg));
                proc.arbitrary_filter_enabled = true;
            } else if (arg == "--binaural") {
                binaural.enabled = true;
                binaural.beat_frequency_hz = std::stof(need_value(i, arg));
            } else if (arg == "--binaural-mono") {
                binaural.mono = std::stof(need_value(i, arg));
            } else if (arg == "--binaural-mode") {
                binaural.stereo_mode = parse_binaural_mode(need_value(i, arg));
            } else if (arg == "--binaural-envelope") {
                binaural_env = parse_breakpoints(need_value(i, arg));
            } else if (arg == "-h" || arg == "--help") {
                print_usage(argv[0]);
                return EXIT_SUCCESS;
            } else if (!arg.empty() && arg[0] == '-') {
                throw std::runtime_error("unknown option '" + arg + "'");
            } else if (input_path.empty()) {
                input_path = arg;
            } else if (output_path.empty()) {
                output_path = arg;
            } else {
                throw std::runtime_error("unexpected argument '" + arg + "'");
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << "\n";
        print_usage(argv[0]);
        return EXIT_FAILURE;
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

    auto num_channels_in = wav_in.channels;
    auto sample_rate = wav_in.sampleRate;
    auto total_frames = wav_in.totalPCMFrameCount;

    std::vector<float> interleaved(total_frames * num_channels_in);
    drwav_read_pcm_frames_f32(&wav_in, total_frames, interleaved.data());
    drwav_uninit(&wav_in);

    // Deinterleave
    std::vector<float> left(total_frames), right(total_frames);
    if (num_channels_in == 1) {
        left = std::move(interleaved);
        right = left;
    } else {
        for (drwav_uint64 i = 0; i < total_frames; i++) {
            left[i] = interleaved[i * 2];
            right[i] = interleaved[i * 2 + 1];
        }
    }

    opts.sample_rate = static_cast<float>(sample_rate);

    std::cerr << "Input: " << total_frames << " frames, " << num_channels_in
              << "ch, " << sample_rate << " Hz\n"
              << "Stretch: " << opts.stretch << "x, FFT: " << opts.fft_size
              << ", Backend: " << paulstretch::fft_backend_name() << "\n";

    paulstretch::OfflineRenderer renderer(opts);
    renderer.set_process_options(proc);
    if (!stretch_env.empty()) renderer.set_stretch_envelope(stretch_env);
    if (!arb_filter.empty()) renderer.set_arbitrary_filter(arb_filter);

    bool output_stereo = (num_channels_in == 2) || binaural.enabled;

    std::vector<float> out_left, out_right;
    if (num_channels_in == 2) {
        auto rendered = renderer.render_stereo(left, right);
        out_left = std::move(rendered.left);
        out_right = std::move(rendered.right);
    } else {
        out_left = renderer.render_mono(left);
        if (output_stereo) out_right = out_left;
    }

    auto out_frames = out_left.size();
    std::cerr << "Output: " << out_frames << " frames\n";

    // Optional binaural-beats post-process
    if (binaural.enabled) {
        paulstretch::BinauralBeatsProcessor bp(static_cast<float>(sample_rate));
        bp.set_options(binaural);
        if (!binaural_env.empty()) bp.set_frequency_envelope(binaural_env);

        // Process in chunks so the frequency envelope can evolve along the output.
        constexpr int kChunk = 4096;
        for (size_t pos = 0; pos < out_frames; pos += kChunk) {
            int n = static_cast<int>(std::min<size_t>(kChunk, out_frames - pos));
            float position_pct = out_frames > 0
                ? 100.0f * static_cast<float>(pos) / static_cast<float>(out_frames)
                : 0.0f;
            bp.process(out_left.data() + pos, out_right.data() + pos, n, position_pct);
        }
        std::cerr << "Binaural: " << binaural.beat_frequency_hz << " Hz\n";
    }

    // Write output
    drwav wav_out;
    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = output_stereo ? 2 : 1;
    format.sampleRate = sample_rate;
    format.bitsPerSample = 32;

    if (!drwav_init_file_write(&wav_out, output_path.c_str(), &format, nullptr)) {
        std::cerr << "error: cannot create '" << output_path << "'\n";
        return EXIT_FAILURE;
    }

    if (output_stereo) {
        std::vector<float> out_interleaved(out_frames * 2);
        for (size_t i = 0; i < out_frames; i++) {
            out_interleaved[i * 2] = out_left[i];
            out_interleaved[i * 2 + 1] = out_right[i];
        }
        drwav_write_pcm_frames(&wav_out, out_frames, out_interleaved.data());
    } else {
        drwav_write_pcm_frames(&wav_out, out_frames, out_left.data());
    }
    drwav_uninit(&wav_out);

    std::cerr << "Wrote " << output_path << "\n";
    return EXIT_SUCCESS;
}
