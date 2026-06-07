#include "paulstretch/paulstretch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>

#if defined(PAULSTRETCH_USE_PFFFT)
#include <pffft/pffft.h>
#elif defined(PAULSTRETCH_USE_ACCELERATE)
#include <Accelerate/Accelerate.h>
#elif defined(KISSFFT)
#include <kiss_fftr.h>
#else
#include <fftw3.h>
#endif

namespace paulstretch {
namespace {

constexpr float kPi = 3.14159265358979323846f;

// ── Envelope interpolation ─────────────────────────────────────────────────

float envelope_at(const std::vector<Breakpoint> &env, float position) {
    if (env.empty()) return 1.0f;
    if (position <= env.front().position) return env.front().value;
    if (position >= env.back().position) return env.back().value;
    auto it = std::lower_bound(env.begin(), env.end(), position,
        [](const Breakpoint &bp, float pos) { return bp.position < pos; });
    if (it == env.begin()) return it->value;
    auto prev = std::prev(it);
    float t = (position - prev->position) / (it->position - prev->position);
    return prev->value + t * (it->value - prev->value);
}

void sort_breakpoints(std::vector<Breakpoint> &points) {
    std::sort(points.begin(), points.end(),
              [](const Breakpoint &a, const Breakpoint &b) { return a.position < b.position; });
}

// ── FFT ────────────────────────────────────────────────────────────────────

class FFT {
public:
    float *smp;
    float *freq;
    int nsamples;

    explicit FFT(int nsamples_) : nsamples(nsamples_) {
        if (nsamples % 2 != 0) ++nsamples;
#if defined(PAULSTRETCH_USE_ACCELERATE)
        { int p = 1; while (p < nsamples) p <<= 1; nsamples = p; }
#endif

        smp_storage_.assign(nsamples, 0.0f);
        freq_storage_.assign(nsamples / 2 + 1, 0.0f);
        window_storage_.assign(nsamples, 0.707f);
        smp = smp_storage_.data();
        freq = freq_storage_.data();
        window_.data = window_storage_.data();
        window_.type = Window::Rectangular;

#if defined(PAULSTRETCH_USE_PFFFT)
        if (!pffft_is_valid_size(nsamples, PFFFT_REAL))
            throw std::invalid_argument("invalid PFFFT transform size");
        pffft_setup_.reset(pffft_new_setup(nsamples, PFFFT_REAL));
        if (!pffft_setup_) throw std::runtime_error("PFFFT setup failed");
        pffft_input_.reset(static_cast<float *>(pffft_aligned_malloc(sizeof(float) * nsamples)));
        pffft_output_.reset(static_cast<float *>(pffft_aligned_malloc(sizeof(float) * nsamples)));
        pffft_work_.reset(static_cast<float *>(pffft_aligned_malloc(sizeof(float) * nsamples)));
        if (!pffft_input_ || !pffft_output_ || !pffft_work_) throw std::bad_alloc();
        std::fill_n(pffft_input_.get(), nsamples, 0.0f);
        std::fill_n(pffft_output_.get(), nsamples, 0.0f);
        std::fill_n(pffft_work_.get(), nsamples, 0.0f);
#elif defined(PAULSTRETCH_USE_ACCELERATE)
        log2n_ = static_cast<vDSP_Length>(std::round(std::log2(nsamples)));
        fft_setup_ = vDSP_create_fftsetup(log2n_, FFT_RADIX2);
        if (!fft_setup_) throw std::runtime_error("vDSP_create_fftsetup failed");
        split_realp_.assign(nsamples / 2, 0.0f);
        split_imagp_.assign(nsamples / 2, 0.0f);
        split_.realp = split_realp_.data();
        split_.imagp = split_imagp_.data();
#elif defined(KISSFFT)
        real_buf_.assign(nsamples + 2, 0.0f);
        cpx_buf_.resize(nsamples / 2 + 2);
        for (auto &c : cpx_buf_) c.r = c.i = 0.0f;
        fwd_cfg_ = kiss_fftr_alloc(nsamples, 0, 0, 0);
        inv_cfg_ = kiss_fftr_alloc(nsamples, 1, 0, 0);
#else
        data_buf_.assign(nsamples, 0.0f);
        fwd_plan_ = fftwf_plan_r2r_1d(nsamples, data_buf_.data(), data_buf_.data(), FFTW_R2HC, FFTW_ESTIMATE);
        inv_plan_ = fftwf_plan_r2r_1d(nsamples, data_buf_.data(), data_buf_.data(), FFTW_HC2R, FFTW_ESTIMATE);
#endif
        rand_seed_ = start_rand_seed_;
        start_rand_seed_ += 161103;
    }

    ~FFT() {
#if defined(PAULSTRETCH_USE_PFFFT)
        // unique_ptrs handle cleanup
#elif defined(PAULSTRETCH_USE_ACCELERATE)
        vDSP_destroy_fftsetup(fft_setup_);
#elif defined(KISSFFT)
        free(fwd_cfg_);
        free(inv_cfg_);
#else
        fftwf_destroy_plan(fwd_plan_);
        fftwf_destroy_plan(inv_plan_);
#endif
    }

    FFT(const FFT &) = delete;
    FFT &operator=(const FFT &) = delete;

    void smp2freq() {
#if defined(PAULSTRETCH_USE_PFFFT)
        std::copy_n(smp, nsamples, pffft_input_.get());
        pffft_transform_ordered(pffft_setup_.get(), pffft_input_.get(),
                                pffft_output_.get(), pffft_work_.get(), PFFFT_FORWARD);
#elif defined(PAULSTRETCH_USE_ACCELERATE)
        vDSP_ctoz(reinterpret_cast<const DSPComplex *>(smp), 2, &split_, 1,
                  static_cast<vDSP_Length>(nsamples / 2));
        vDSP_fft_zrip(fft_setup_, &split_, 1, log2n_, FFT_FORWARD);
#elif defined(KISSFFT)
        std::copy_n(smp, nsamples, real_buf_.data());
        kiss_fftr(fwd_cfg_, real_buf_.data(), cpx_buf_.data());
#else
        std::copy_n(smp, nsamples, data_buf_.data());
        fftwf_execute(fwd_plan_);
#endif
        for (int i = 1; i < nsamples / 2; i++) {
#if defined(PAULSTRETCH_USE_PFFFT)
            float c = pffft_output_.get()[2 * i];
            float s = pffft_output_.get()[2 * i + 1];
#elif defined(PAULSTRETCH_USE_ACCELERATE)
            float c = split_.realp[i];
            float s = split_.imagp[i];
#elif defined(KISSFFT)
            float c = cpx_buf_[i].r;
            float s = cpx_buf_[i].i;
#else
            float c = data_buf_[i];
            float s = data_buf_[nsamples - i];
#endif
            freq[i] = std::sqrt(c * c + s * s);
        }
        freq[0] = 0.0f;
        freq[nsamples / 2] = 0.0f;
    }

    void freq2smp() {
        float inv_2p15_2pi = 1.0f / 16384.0f * kPi;
#if defined(PAULSTRETCH_USE_PFFFT)
        std::fill_n(pffft_input_.get(), nsamples, 0.0f);
#elif defined(PAULSTRETCH_USE_ACCELERATE)
        std::fill_n(split_.realp, nsamples / 2, 0.0f);
        std::fill_n(split_.imagp, nsamples / 2, 0.0f);
#endif
        for (int i = 1; i < nsamples / 2; i++) {
            rand_seed_ = rand_seed_ * 1103515245 + 12345;
            unsigned int r = (rand_seed_ >> 16) & 0x7fff;
            float phase = r * inv_2p15_2pi;
#if defined(PAULSTRETCH_USE_PFFFT)
            pffft_input_.get()[2 * i] = freq[i] * std::cos(phase);
            pffft_input_.get()[2 * i + 1] = freq[i] * std::sin(phase);
#elif defined(PAULSTRETCH_USE_ACCELERATE)
            split_.realp[i] = freq[i] * std::cos(phase);
            split_.imagp[i] = freq[i] * std::sin(phase);
#elif defined(KISSFFT)
            cpx_buf_[i].r = freq[i] * std::cos(phase);
            cpx_buf_[i].i = freq[i] * std::sin(phase);
#else
            data_buf_[i] = freq[i] * std::cos(phase);
            data_buf_[nsamples - i] = freq[i] * std::sin(phase);
#endif
        }

#if defined(PAULSTRETCH_USE_PFFFT)
        pffft_transform_ordered(pffft_setup_.get(), pffft_input_.get(),
                                pffft_output_.get(), pffft_work_.get(), PFFFT_BACKWARD);
        for (int i = 0; i < nsamples; i++) smp[i] = pffft_output_.get()[i] / nsamples;
#elif defined(PAULSTRETCH_USE_ACCELERATE)
        split_.realp[0] = 0.0f;
        split_.imagp[0] = 0.0f;
        vDSP_fft_zrip(fft_setup_, &split_, 1, log2n_, FFT_INVERSE);
        vDSP_ztoc(&split_, 1, reinterpret_cast<DSPComplex *>(smp), 2,
                  static_cast<vDSP_Length>(nsamples / 2));
        const float scale = 1.0f / static_cast<float>(2 * nsamples);
        vDSP_vsmul(smp, 1, &scale, smp, 1, static_cast<vDSP_Length>(nsamples));
#elif defined(KISSFFT)
        cpx_buf_[0].r = cpx_buf_[0].i = 0.0f;
        cpx_buf_[nsamples / 2].r = cpx_buf_[nsamples / 2].i = 0.0f;
        kiss_fftri(inv_cfg_, cpx_buf_.data(), real_buf_.data());
        for (int i = 0; i < nsamples; i++) smp[i] = real_buf_[i] / nsamples;
#else
        data_buf_[0] = data_buf_[nsamples / 2] = 0.0f;
        if (nsamples / 2 + 1 < nsamples) data_buf_[nsamples / 2 + 1] = 0.0f;
        fftwf_execute(inv_plan_);
        for (int i = 0; i < nsamples; i++) smp[i] = data_buf_[i] / nsamples;
#endif
    }

    void applywindow(Window type) {
        if (window_.type != type) {
            window_.type = type;
            for (int i = 0; i < nsamples; i++) {
                switch (type) {
                case Window::Rectangular:
                    window_.data[i] = 0.707f;
                    break;
                case Window::Hamming:
                    window_.data[i] = 0.53836f - 0.46164f * std::cos(2 * kPi * i / (nsamples + 1.0f));
                    break;
                case Window::Hann:
                    window_.data[i] = 0.5f * (1.0f - std::cos(2 * kPi * i / (nsamples - 1.0f)));
                    break;
                case Window::Blackman:
                    window_.data[i] = 0.42f - 0.5f * std::cos(2 * kPi * i / (nsamples - 1.0f))
                                     + 0.08f * std::cos(4 * kPi * i / (nsamples - 1.0f));
                    break;
                case Window::BlackmanHarris:
                    window_.data[i] = 0.35875f - 0.48829f * std::cos(2 * kPi * i / (nsamples - 1.0f))
                                     + 0.14128f * std::cos(4 * kPi * i / (nsamples - 1.0f))
                                     - 0.01168f * std::cos(6 * kPi * i / (nsamples - 1.0f));
                    break;
                }
            }
        }
        for (int i = 0; i < nsamples; i++) smp[i] *= window_.data[i];
    }

private:
    struct { float *data; Window type; } window_;
    std::vector<float> smp_storage_;
    std::vector<float> freq_storage_;
    std::vector<float> window_storage_;
    unsigned int rand_seed_;
    static inline unsigned int start_rand_seed_ = 1;

#if defined(PAULSTRETCH_USE_PFFFT)
    struct PffftSetupDel { void operator()(PFFFT_Setup *s) const noexcept { if (s) pffft_destroy_setup(s); } };
    struct PffftFreeDel { void operator()(float *p) const noexcept { if (p) pffft_aligned_free(p); } };
    std::unique_ptr<PFFFT_Setup, PffftSetupDel> pffft_setup_;
    std::unique_ptr<float, PffftFreeDel> pffft_input_;
    std::unique_ptr<float, PffftFreeDel> pffft_output_;
    std::unique_ptr<float, PffftFreeDel> pffft_work_;
#elif defined(PAULSTRETCH_USE_ACCELERATE)
    FFTSetup fft_setup_{};
    vDSP_Length log2n_{};
    DSPSplitComplex split_{};
    std::vector<float> split_realp_;
    std::vector<float> split_imagp_;
#elif defined(KISSFFT)
    kiss_fftr_cfg fwd_cfg_, inv_cfg_;
    std::vector<kiss_fft_scalar> real_buf_;
    std::vector<kiss_fft_cpx> cpx_buf_;
#else
    fftwf_plan fwd_plan_, inv_plan_;
    std::vector<float> data_buf_;
#endif
};

// ── Original Win32 spectral post-processing chain ──────────────────────────

class SpectralProcessor {
public:
    void set_options(ProcessOptions options) { options_ = options; }
    const ProcessOptions &options() const { return options_; }

    void set_arbitrary_filter(std::vector<Breakpoint> filter) {
        sort_breakpoints(filter);
        arbitrary_filter_ = std::move(filter);
    }

    void clear_arbitrary_filter() { arbitrary_filter_.clear(); }
    const std::vector<Breakpoint> &arbitrary_filter() const { return arbitrary_filter_; }

    void process(float *freq, int nfreq, float sample_rate) {
        if (!has_active_options()) return;
        ensure_size(nfreq);

        if (options_.harmonics_enabled) {
            copy(freq, input_.data(), nfreq);
            do_harmonics(input_.data(), freq, nfreq, sample_rate);
        }
        if (options_.tonal_noise_enabled) {
            copy(freq, input_.data(), nfreq);
            do_tonal_vs_noise(input_.data(), freq, nfreq, sample_rate);
        }
        if (options_.frequency_shift_enabled) {
            copy(freq, input_.data(), nfreq);
            do_freq_shift(input_.data(), freq, nfreq, sample_rate);
        }
        if (options_.pitch_shift_enabled) {
            copy(freq, input_.data(), nfreq);
            do_pitch_shift(input_.data(), freq, nfreq,
                           std::pow(2.0f, options_.pitch_shift_cents / 1200.0f));
        }
        if (options_.octave_enabled) {
            copy(freq, input_.data(), nfreq);
            do_octave(input_.data(), freq, nfreq);
        }
        if (options_.spread_enabled) {
            copy(freq, input_.data(), nfreq);
            do_spread(input_.data(), freq, nfreq, sample_rate, options_.spread_bandwidth);
        }
        if (options_.filter_enabled) {
            copy(freq, input_.data(), nfreq);
            do_filter(input_.data(), freq, nfreq, sample_rate);
        }
        if (options_.arbitrary_filter_enabled && !arbitrary_filter_.empty()) {
            copy(freq, input_.data(), nfreq);
            do_arbitrary_filter(input_.data(), freq, nfreq, sample_rate);
        }
        if (options_.compressor_enabled) {
            copy(freq, input_.data(), nfreq);
            do_compressor(input_.data(), freq, nfreq);
        }
    }

private:
    bool has_active_options() const {
        return options_.harmonics_enabled ||
               options_.tonal_noise_enabled ||
               options_.frequency_shift_enabled ||
               options_.pitch_shift_enabled ||
               options_.octave_enabled ||
               options_.spread_enabled ||
               options_.filter_enabled ||
               (options_.arbitrary_filter_enabled && !arbitrary_filter_.empty()) ||
               options_.compressor_enabled;
    }

    void ensure_size(int nfreq) {
        if (static_cast<int>(input_.size()) == nfreq) return;
        input_.assign(nfreq, 0.0f);
        sum_.assign(nfreq, 0.0f);
        tmp1_.assign(nfreq, 0.0f);
        tmp2_.assign(nfreq, 0.0f);
    }

    static void copy(const float *src, float *dst, int nfreq) {
        std::copy_n(src, nfreq, dst);
    }

    static void add(float *dst, const float *src, float scale, int nfreq) {
        for (int i = 0; i < nfreq; i++) dst[i] += src[i] * scale;
    }

    static void zero(float *dst, int nfreq) {
        std::fill_n(dst, nfreq, 0.0f);
    }

    static float profile(float fi, float bwi) {
        const float x = (fi / bwi) * (fi / bwi);
        if (x > 14.71280603f) return 0.0f;
        return std::exp(-x);
    }

    void do_harmonics(const float *freq1, float *freq2, int nfreq, float sample_rate) {
        float fundamental = options_.harmonics_frequency_hz;
        const float bandwidth = options_.harmonics_bandwidth_cents;
        const int harmonics = std::max(1, options_.harmonics_count);
        if (fundamental < 10.0f) fundamental = 10.0f;

        float *amp = tmp1_.data();
        zero(amp, nfreq);

        for (int nh = 1; nh <= harmonics; nh++) {
            const float harmonic_hz = nh * fundamental;
            if (harmonic_hz >= sample_rate * 0.5f) break;

            const float bw_hz = (std::pow(2.0f, bandwidth / 1200.0f) - 1.0f) * harmonic_hz;
            const float bwi = bw_hz / (2.0f * sample_rate);
            if (bwi <= 0.0f) continue;
            const float fi = harmonic_hz / sample_rate;

            for (int i = 1; i < nfreq; i++) {
                amp[i] += profile((i / static_cast<float>(nfreq) * 0.5f) - fi, bwi);
            }
        }

        float max_amp = 0.0f;
        for (int i = 1; i < nfreq; i++) max_amp = std::max(max_amp, amp[i]);
        if (max_amp < 1e-8f) max_amp = 1e-8f;

        for (int i = 1; i < nfreq; i++) {
            float a = amp[i] / max_amp;
            if (!options_.harmonics_gauss) a = (a < 0.368f ? 0.0f : 1.0f);
            freq2[i] = freq1[i] * a;
        }
        if (nfreq > 0) freq2[0] = 0.0f;
    }

    void do_freq_shift(const float *freq1, float *freq2, int nfreq, float sample_rate) {
        zero(freq2, nfreq);
        const int bin_shift = static_cast<int>(
            options_.frequency_shift_hz / (sample_rate * 0.5f) * nfreq);
        for (int i = 0; i < nfreq; i++) {
            const int dst = bin_shift + i;
            if (dst > 0 && dst < nfreq) freq2[dst] = freq1[i];
        }
    }

    static void do_pitch_shift(const float *freq1, float *freq2, int nfreq, float ratio) {
        zero(freq2, nfreq);
        if (ratio < 1.0f) {
            for (int i = 0; i < nfreq; i++) {
                const int dst = static_cast<int>(i * ratio);
                if (dst >= nfreq) break;
                freq2[dst] += freq1[i];
            }
        } else {
            const float inv_ratio = 1.0f / ratio;
            for (int i = 0; i < nfreq; i++) {
                const int src = std::min(nfreq - 1, static_cast<int>(i * inv_ratio));
                freq2[i] = freq1[src];
            }
        }
    }

    void do_octave(const float *freq1, float *freq2, int nfreq) {
        zero(sum_.data(), nfreq);
        if (options_.octave_minus2 > 1e-3f) {
            do_pitch_shift(freq1, tmp1_.data(), nfreq, 0.25f);
            add(sum_.data(), tmp1_.data(), options_.octave_minus2, nfreq);
        }
        if (options_.octave_minus1 > 1e-3f) {
            do_pitch_shift(freq1, tmp1_.data(), nfreq, 0.5f);
            add(sum_.data(), tmp1_.data(), options_.octave_minus1, nfreq);
        }
        if (options_.octave_0 > 1e-3f) {
            add(sum_.data(), freq1, options_.octave_0, nfreq);
        }
        if (options_.octave_plus1 > 1e-3f) {
            do_pitch_shift(freq1, tmp1_.data(), nfreq, 2.0f);
            add(sum_.data(), tmp1_.data(), options_.octave_plus1, nfreq);
        }
        if (options_.octave_plus15 > 1e-3f) {
            do_pitch_shift(freq1, tmp1_.data(), nfreq, 3.0f);
            add(sum_.data(), tmp1_.data(), options_.octave_plus15, nfreq);
        }
        if (options_.octave_plus2 > 1e-3f) {
            do_pitch_shift(freq1, tmp1_.data(), nfreq, 4.0f);
            add(sum_.data(), tmp1_.data(), options_.octave_plus2, nfreq);
        }

        float sum = 0.01f + options_.octave_minus2 + options_.octave_minus1 +
                    options_.octave_0 + options_.octave_plus1 +
                    options_.octave_plus15 + options_.octave_plus2;
        if (sum < 0.5f) sum = 0.5f;
        for (int i = 0; i < nfreq; i++) freq2[i] = sum_[i] / sum;
    }

    void do_filter(const float *freq1, float *freq2, int nfreq, float sample_rate) {
        const float low = std::min(options_.filter_low_hz, options_.filter_high_hz);
        const float high = std::max(options_.filter_low_hz, options_.filter_high_hz);
        const int ilow = static_cast<int>(low / sample_rate * nfreq * 2.0f);
        const int ihigh = static_cast<int>(high / sample_rate * nfreq * 2.0f);
        float damp = 1.0f;
        const float damp_ratio = 1.0f - std::pow(options_.filter_high_damp * 0.5f, 4.0f);
        for (int i = 0; i < nfreq; i++) {
            float a = (i >= ilow && i < ihigh) ? 1.0f : 0.0f;
            if (options_.filter_stop) a = 1.0f - a;
            freq2[i] = freq1[i] * a * damp;
            damp *= damp_ratio + 1e-8f;
        }
    }

    void do_arbitrary_filter(const float *freq1, float *freq2, int nfreq, float sample_rate) {
        constexpr float min_freq = 20.0f;
        constexpr float max_freq = 25000.0f;
        const float log_ratio = std::log(max_freq / min_freq);
        for (int i = 0; i < nfreq; i++) {
            const float hz = i / static_cast<float>(nfreq) * sample_rate * 0.5f;
            float x = 0.0f;
            if (hz > min_freq && log_ratio > 0.0f)
                x = std::log(hz / min_freq) / log_ratio;
            x = std::clamp(x, 0.0f, 1.0f);
            freq2[i] = freq1[i] * envelope_at(arbitrary_filter_, x);
        }
    }

    void do_spread(const float *freq1, float *freq2, int nfreq, float sample_rate, float bandwidth) {
        const float min_freq = 20.0f;
        const float max_freq = 0.5f * sample_rate;
        if (max_freq <= min_freq || nfreq <= 1) {
            copy(freq1, freq2, nfreq);
            return;
        }

        const float log_minfreq = std::log(min_freq);
        const float log_maxfreq = std::log(max_freq);
        for (int i = 0; i < nfreq; i++) {
            const float freqx = i / static_cast<float>(nfreq);
            const float x = std::exp(log_minfreq + freqx * (log_maxfreq - log_minfreq))
                          / max_freq * nfreq;
            float y = 0.0f;
            if (x < nfreq) {
                const int x0 = std::min(nfreq - 1, static_cast<int>(std::floor(x)));
                const int x1 = std::min(nfreq - 1, x0 + 1);
                const float xp = x - x0;
                y = freq1[x0] * (1.0f - xp) + freq1[x1] * xp;
            }
            tmp1_[i] = y;
        }

        const int passes = 2;
        float a = 1.0f - std::pow(2.0f, -bandwidth * bandwidth * 10.0f);
        a = std::pow(a, 8192.0f / nfreq * passes);
        for (int k = 0; k < passes; k++) {
            tmp1_[0] = 0.0f;
            for (int i = 1; i < nfreq; i++) tmp1_[i] = tmp1_[i - 1] * a + tmp1_[i] * (1.0f - a);
            tmp1_[nfreq - 1] = 0.0f;
            for (int i = nfreq - 2; i > 0; i--) tmp1_[i] = tmp1_[i + 1] * a + tmp1_[i] * (1.0f - a);
        }

        freq2[0] = 0.0f;
        const float log_maxfreq_d_minfreq = std::log(max_freq / min_freq);
        for (int i = 1; i < nfreq; i++) {
            const float freqx = i / static_cast<float>(nfreq);
            const float x = std::log((freqx * max_freq) / min_freq) / log_maxfreq_d_minfreq * nfreq;
            float y = 0.0f;
            if (x > 0.0f && x < nfreq) {
                const int x0 = std::min(nfreq - 1, static_cast<int>(std::floor(x)));
                const int x1 = std::min(nfreq - 1, x0 + 1);
                const float xp = x - x0;
                y = tmp1_[x0] * (1.0f - xp) + tmp1_[x1] * xp;
            }
            freq2[i] = y;
        }
    }

    void do_compressor(const float *freq1, float *freq2, int nfreq) {
        float rms = 0.0f;
        for (int i = 0; i < nfreq; i++) rms += freq1[i] * freq1[i];
        rms = std::sqrt(rms / nfreq) * 0.1f;
        if (rms < 1e-3f) rms = 1e-3f;

        const float ratio = std::pow(rms, -options_.compressor_power);
        for (int i = 0; i < nfreq; i++) freq2[i] = freq1[i] * ratio;
    }

    void do_tonal_vs_noise(const float *freq1, float *freq2, int nfreq, float sample_rate) {
        do_spread(freq1, tmp1_.data(), nfreq, sample_rate, options_.tonal_noise_bandwidth);

        if (options_.tonal_noise_preserve >= 0.0f) {
            const float mul = std::pow(10.0f, options_.tonal_noise_preserve) - 1.0f;
            for (int i = 0; i < nfreq; i++) {
                const float smooth_x = tmp1_[i] + 1e-6f;
                freq2[i] = std::max(0.0f, freq1[i] - smooth_x * mul);
            }
        } else {
            const float mul = std::pow(5.0f, 1.0f + options_.tonal_noise_preserve) - 1.0f;
            for (int i = 0; i < nfreq; i++) {
                const float smooth_x = tmp1_[i] + 1e-6f;
                const float result = freq1[i] - smooth_x * mul + 0.1f * mul;
                freq2[i] = result < 0.0f ? freq1[i] : 0.0f;
            }
        }
    }

    ProcessOptions options_;
    std::vector<Breakpoint> arbitrary_filter_;
    std::vector<float> input_;
    std::vector<float> sum_;
    std::vector<float> tmp1_;
    std::vector<float> tmp2_;
};

// ── Stretcher ──────────────────────────────────────────────────────────────

class Stretcher {
public:
    float *out_buf;

    Stretcher(float stretch_factor, int fft_size, Window window,
              float sample_rate, int stereo_mode,
              const std::vector<Breakpoint> *envelope)
        : envelope_(envelope), window_(window), stereo_mode_(stereo_mode),
          stretch_factor_(stretch_factor), sample_rate_(sample_rate),
          onset_sensitivity_(0.0f), remained_samples_(0.0),
          extra_onset_time_credit_(0.0), c_pos_percents_(0.0f),
          skip_samples_(0), require_new_buffer_(false), freezing_(false) {
        bufsize_ = fft_size;
        if (bufsize_ < 8) bufsize_ = 8;
#if defined(PAULSTRETCH_USE_PFFFT)
        bufsize_ = pffft_nearest_transform_size(bufsize_ * 2, PFFFT_REAL, 1) / 2;
#elif defined(PAULSTRETCH_USE_ACCELERATE)
        { int p = 1; while (p < bufsize_ * 2) p <<= 1; bufsize_ = p / 2; }
#endif
        out_buf_storage_.assign(bufsize_, 0.0f);
        out_buf = out_buf_storage_.data();
        old_freq_.assign(bufsize_, 0.0f);
        very_old_smps_.assign(bufsize_, 0.0f);
        new_smps_.assign(bufsize_, 0.0f);
        old_smps_.assign(bufsize_, 0.0f);
        old_out_smps_.assign(bufsize_ * 2, 0.0f);

        infft_ = std::make_unique<FFT>(bufsize_ * 2);
        fft_ = std::make_unique<FFT>(bufsize_ * 2);
        outfft_ = std::make_unique<FFT>(bufsize_ * 2);
    }

    ~Stretcher() = default;
    Stretcher(const Stretcher &) = delete;
    Stretcher &operator=(const Stretcher &) = delete;

    int bufsize() const { return bufsize_; }
    int max_bufsize() const { return bufsize_ * 3; }

    void set_onset_detection_sensitivity(float s) {
        onset_sensitivity_ = s;
        if (s < 1e-3f) extra_onset_time_credit_ = 0.0;
    }

    // Hot-swap the stretch factor mid-stream. The next process() picks it up
    // when computing how much input to advance over. No FFT state is reset,
    // so output is continuous across the change.
    void set_stretch_factor(float s) { stretch_factor_ = s; }

    // Hot-swap the envelope pointer. `env` may be nullptr to disable the
    // envelope; otherwise must outlive the next process() call. Updating
    // the storage in-place under the existing pointer also works.
    void set_envelope(const std::vector<Breakpoint> *env) { envelope_ = env; }
    void set_process_options(ProcessOptions options) { spectral_.set_options(options); }
    const ProcessOptions &process_options() const { return spectral_.options(); }
    void set_arbitrary_filter(std::vector<Breakpoint> filter) {
        spectral_.set_arbitrary_filter(std::move(filter));
    }
    void clear_arbitrary_filter() { spectral_.clear_arbitrary_filter(); }
    const std::vector<Breakpoint> &arbitrary_filter() const {
        return spectral_.arbitrary_filter();
    }

    int get_nsamples(float current_pos_percents) {
        if (freezing_) return 0;
        c_pos_percents_ = current_pos_percents;
        return require_new_buffer_ ? bufsize_ : 0;
    }

    // Pure query — returns the same value get_nsamples() would, without
    // mutating the recorded position. For the streaming API where we want
    // to ask "how much input next?" separately from "advance with this
    // position".
    int peek_nsamples() const {
        if (freezing_) return 0;
        return require_new_buffer_ ? bufsize_ : 0;
    }

    int get_nsamples_for_fill() { return max_bufsize(); }

    int get_skip_nsamples() {
        if (freezing_) return 0;
        return skip_samples_;
    }

    void here_is_onset(float onset) {
        if (freezing_) return;
        if (onset > 0.5f) {
            require_new_buffer_ = true;
            extra_onset_time_credit_ += 1.0 - remained_samples_;
            remained_samples_ = 0.0;
            skip_samples_ = 0;
        }
    }

    float process(const float *smps, int nsmps) {
        float onset = 0.0f;

        if (smps != nullptr) {
            if (nsmps != 0 && nsmps != bufsize_ && nsmps != max_bufsize()) return 0.0f;

            if (nsmps != 0) {
                do_analyse_inbuf(smps);
                if (nsmps == max_bufsize()) {
                    for (int k = bufsize_; k < max_bufsize(); k += bufsize_)
                        do_analyse_inbuf(smps + k);
                }
                if (onset_sensitivity_ > 1e-3f) onset = do_detect_onset();
            }

            if (nsmps != 0) {
                do_next_inbuf_smps(smps);
                if (nsmps == max_bufsize()) {
                    for (int k = bufsize_; k < max_bufsize(); k += bufsize_)
                        do_next_inbuf_smps(smps + k);
                }
            }

            // construct the input FFT
            int start_pos = (int)(std::floor(remained_samples_ * bufsize_));
            if (start_pos >= bufsize_) start_pos = bufsize_ - 1;
            for (int i = 0; i < bufsize_ - start_pos; i++)
                fft_->smp[i] = very_old_smps_[i + start_pos];
            for (int i = 0; i < bufsize_; i++)
                fft_->smp[i + bufsize_ - start_pos] = old_smps_[i];
            for (int i = 0; i < start_pos; i++)
                fft_->smp[i + 2 * bufsize_ - start_pos] = new_smps_[i];

            fft_->applywindow(window_);
            fft_->smp2freq();
            for (int i = 0; i < bufsize_; i++) outfft_->freq[i] = fft_->freq[i];
            spectral_.process(outfft_->freq, bufsize_, sample_rate_);

            outfft_->freq2smp();

            // overlap-add with amplitude correction
            float tmp = 1.0f / (float)bufsize_ * kPi;
            float hinv_sqrt2 = 0.853553390593f;
            float ampfactor = 2.0f;

            for (int i = 0; i < bufsize_; i++) {
                float a = 0.5f + 0.5f * std::cos(i * tmp);
                float out = outfft_->smp[i + bufsize_] * (1.0f - a) + old_out_smps_[i] * a;
                out_buf[i] = out * (hinv_sqrt2 - (1.0f - hinv_sqrt2) * std::cos(i * 2.0f * tmp)) * ampfactor;
            }

            for (int i = 0; i < bufsize_ * 2; i++) old_out_smps_[i] = outfft_->smp[i];
        }

        if (!freezing_) {
            long double used_rap = stretch_factor_ * get_stretch_multiplier(c_pos_percents_);
            long double r = 1.0 / used_rap;
            if (extra_onset_time_credit_ > 0) {
                float credit_get = 0.5f * (float)r;
                extra_onset_time_credit_ -= credit_get;
                if (extra_onset_time_credit_ < 0.0) extra_onset_time_credit_ = 0.0;
                r -= credit_get;
            }
            remained_samples_ += r;
            if (remained_samples_ >= 1.0) {
                skip_samples_ = (int)(std::floor(remained_samples_ - 1.0) * bufsize_);
                remained_samples_ = remained_samples_ - std::floor(remained_samples_);
                require_new_buffer_ = true;
            } else {
                require_new_buffer_ = false;
            }
        }

        return onset;
    }

private:
    float get_stretch_multiplier(float pos_percents) {
        if (!envelope_ || envelope_->empty()) return 1.0f;
        return envelope_at(*envelope_, pos_percents / 100.0f);
    }

    void do_analyse_inbuf(const float *smps) {
        for (int i = 0; i < bufsize_; i++) {
            infft_->smp[i] = old_smps_[i];
            infft_->smp[i + bufsize_] = smps[i];
            old_freq_[i] = infft_->freq[i];
        }
        infft_->applywindow(window_);
        infft_->smp2freq();
    }

    void do_next_inbuf_smps(const float *smps) {
        for (int i = 0; i < bufsize_; i++) {
            very_old_smps_[i] = old_smps_[i];
            old_smps_[i] = new_smps_[i];
            new_smps_[i] = smps[i];
        }
    }

    float do_detect_onset() {
        float result = 0.0f;
        if (onset_sensitivity_ > 1e-3f) {
            float os = 0.0f, osinc = 0.0f;
            float osincold = 1e-5f;
            int maxk = 1 + (int)(bufsize_ * 500.0f / (sample_rate_ * 0.5f));
            int k = 0;
            for (int i = 0; i < bufsize_; i++) {
                osinc += infft_->freq[i] - old_freq_[i];
                osincold += old_freq_[i];
                if (k >= maxk) {
                    k = 0;
                    os += osinc / osincold;
                    osinc = 0;
                }
                k++;
            }
            os += osinc;
            if (os < 0.0f) os = 0.0f;

            float os_strength = std::pow(20.0f, 1.0f - onset_sensitivity_) - 1.0f;
            float os_strength_h = os_strength * 0.75f;
            if (os > os_strength_h) {
                result = (os - os_strength_h) / (os_strength - os_strength_h);
                if (result > 1.0f) result = 1.0f;
            }
            if (result > 1.0f) result = 1.0f;
        }
        return result;
    }

    const std::vector<Breakpoint> *envelope_;
    Window window_;
    int stereo_mode_;
    int bufsize_;
    float stretch_factor_;
    float sample_rate_;
    float onset_sensitivity_;

    std::vector<float> out_buf_storage_;
    std::vector<float> old_out_smps_;
    std::vector<float> old_freq_;
    std::vector<float> new_smps_;
    std::vector<float> old_smps_;
    std::vector<float> very_old_smps_;

    std::unique_ptr<FFT> infft_;
    std::unique_ptr<FFT> outfft_;
    std::unique_ptr<FFT> fft_;
    SpectralProcessor spectral_;

    long double remained_samples_;
    long double extra_onset_time_credit_;
    float c_pos_percents_;
    int skip_samples_;
    bool require_new_buffer_;
    bool freezing_;
};

// ── Rendering helpers ──────────────────────────────────────────────────────

RenderOptions sanitize_options(RenderOptions options) {
    if (options.stretch <= 0.0f) throw std::invalid_argument("stretch must be positive");
    if (options.fft_size < 8) options.fft_size = 8;
    if (options.fft_size % 2 != 0) ++options.fft_size;
#if defined(PAULSTRETCH_USE_PFFFT)
    options.fft_size = pffft_nearest_transform_size(options.fft_size * 2, PFFFT_REAL, 1) / 2;
#elif defined(PAULSTRETCH_USE_ACCELERATE)
    { int p = 1; while (p < options.fft_size * 2) p <<= 1; options.fft_size = p / 2; }
#endif
    if (options.sample_rate <= 0.0f) throw std::invalid_argument("sample_rate must be positive");
    options.onset_detection_sensitivity = std::clamp(options.onset_detection_sensitivity, 0.0f, 1.0f);
    return options;
}

std::size_t clamp_advance(std::size_t cursor, int delta, std::size_t limit) {
    if (delta <= 0) return cursor;
    return std::min(cursor + static_cast<std::size_t>(delta), limit);
}

} // anonymous namespace

// ── StreamingStretcher impl ────────────────────────────────────────────────

struct StreamingStretcher::Impl {
    RenderOptions options;
    std::vector<Breakpoint> envelope;
    ProcessOptions process_options;
    std::vector<Breakpoint> arbitrary_filter;
    std::unique_ptr<Stretcher> stretch;
    float zero_input_sample = 0.0f;
    // `first_step` is true until the initial fill step has run. Drives
    // next_input_size()'s "max_input_chunk vs. {0, bufsize}" branch.
    bool first_step = true;
    // Last reported skip_after_step() value, updated each step().
    int skip_after = 0;

    explicit Impl(RenderOptions opts)
        : options(sanitize_options(opts)) {
        rebuild();
    }

    void rebuild() {
        stretch = std::make_unique<Stretcher>(
            options.stretch, options.fft_size, options.window, options.sample_rate,
            /*stereo_mode=*/0,
            envelope.empty() ? nullptr : &envelope);
        stretch->set_onset_detection_sensitivity(options.onset_detection_sensitivity);
        stretch->set_process_options(process_options);
        stretch->set_arbitrary_filter(arbitrary_filter);
        first_step = true;
        skip_after = 0;
    }
};

StreamingStretcher::StreamingStretcher(RenderOptions options)
    : impl_(std::make_unique<Impl>(options)) {}

StreamingStretcher::~StreamingStretcher() = default;

const RenderOptions &StreamingStretcher::options() const { return impl_->options; }

int StreamingStretcher::bufsize() const { return impl_->stretch->bufsize(); }

int StreamingStretcher::max_input_chunk() const { return impl_->stretch->max_bufsize(); }

int StreamingStretcher::next_input_size() const {
    if (impl_->first_step) return impl_->stretch->max_bufsize();
    return impl_->stretch->peek_nsamples();
}

int StreamingStretcher::skip_after_step() const { return impl_->skip_after; }

float StreamingStretcher::step(const float *input, float position_pct, float *output) {
    const float onset = step_without_onset_feedback(input, position_pct, output);
    apply_onset(onset);
    return onset;
}

float StreamingStretcher::step_without_onset_feedback(
    const float *input, float position_pct, float *output) {
    Stretcher &s = *impl_->stretch;
    // Always record position so the envelope is evaluated at the caller's
    // current cursor — important for seek to land on the right envelope
    // value on the very next step.
    const int natural_n = s.get_nsamples(position_pct);
    const int n = impl_->first_step ? s.max_bufsize() : natural_n;
    impl_->first_step = false;

    if (n > 0 && input == nullptr)
        throw std::invalid_argument("input is required when next_input_size() > 0");

    const float *process_input = n > 0 ? input : &impl_->zero_input_sample;
    const float onset = s.process(process_input, n);
    std::copy_n(s.out_buf, s.bufsize(), output);
    impl_->skip_after = s.get_skip_nsamples();
    return onset;
}

void StreamingStretcher::apply_onset(float onset) {
    impl_->stretch->here_is_onset(onset);
    impl_->skip_after = impl_->stretch->get_skip_nsamples();
}

void StreamingStretcher::set_stretch_envelope(std::vector<Breakpoint> envelope) {
    sort_breakpoints(envelope);
    impl_->envelope = std::move(envelope);
    // Hot-swap: update the inner Stretcher's envelope pointer in place.
    // No DSP state reset, so audio stays continuous across the swap.
    impl_->stretch->set_envelope(impl_->envelope.empty() ? nullptr : &impl_->envelope);
}

void StreamingStretcher::clear_stretch_envelope() {
    impl_->envelope.clear();
    impl_->stretch->set_envelope(nullptr);
}

void StreamingStretcher::set_stretch_factor(float stretch) {
    impl_->options.stretch = stretch;
    impl_->stretch->set_stretch_factor(stretch);
}

const std::vector<Breakpoint> &StreamingStretcher::stretch_envelope() const {
    return impl_->envelope;
}

void StreamingStretcher::set_process_options(ProcessOptions options) {
    impl_->process_options = options;
    impl_->stretch->set_process_options(options);
}

const ProcessOptions &StreamingStretcher::process_options() const {
    return impl_->process_options;
}

void StreamingStretcher::set_arbitrary_filter(std::vector<Breakpoint> filter) {
    sort_breakpoints(filter);
    impl_->arbitrary_filter = std::move(filter);
    impl_->stretch->set_arbitrary_filter(impl_->arbitrary_filter);
}

void StreamingStretcher::clear_arbitrary_filter() {
    impl_->arbitrary_filter.clear();
    impl_->stretch->clear_arbitrary_filter();
}

const std::vector<Breakpoint> &StreamingStretcher::arbitrary_filter() const {
    return impl_->arbitrary_filter;
}

void StreamingStretcher::set_onset_detection_sensitivity(float s) {
    impl_->options.onset_detection_sensitivity = std::clamp(s, 0.0f, 1.0f);
    impl_->stretch->set_onset_detection_sensitivity(impl_->options.onset_detection_sensitivity);
}

void StreamingStretcher::reset() {
    impl_->rebuild();
}

namespace {

// ── Offline render glue (built on StreamingStretcher) ──────────────────────

// Core streaming loop: drives a StreamingStretcher over `input` and delivers
// each `bufsize` output chunk to `sink` (never materialising the full output).
// Both render_channel (whole-buffer) and render_mono_chunked build on this.
void stream_channel(
    const std::vector<float> &input,
    const RenderOptions &options,
    const std::vector<Breakpoint> &envelope,
    const ProcessOptions &process_options,
    const std::vector<Breakpoint> &arbitrary_filter,
    const OfflineRenderer::ChunkSink &sink) {
    if (input.empty()) return;

    StreamingStretcher stretch(options);
    if (!envelope.empty()) stretch.set_stretch_envelope(envelope);
    stretch.set_process_options(process_options);
    if (!arbitrary_filter.empty()) stretch.set_arbitrary_filter(arbitrary_filter);

    const int bufsize = stretch.bufsize();
    std::vector<float> out_chunk(bufsize, 0.0f);
    std::vector<float> in_chunk(stretch.max_input_chunk(), 0.0f);

    std::size_t cursor = 0;
    bool first = true;

    while (true) {
        const float pos_pct = 100.0f * static_cast<float>(cursor) / static_cast<float>(input.size());
        const int want = first ? stretch.max_input_chunk() : stretch.next_input_size();

        if (want > 0 && cursor >= input.size() && !first) break;

        // Gather `want` frames of input (zero-pad past EOF).
        if (want > 0) {
            const std::size_t avail = input.size() - cursor;
            const std::size_t take = std::min<std::size_t>(avail, static_cast<std::size_t>(want));
            std::copy_n(input.data() + cursor, take, in_chunk.data());
            if (take < static_cast<std::size_t>(want))
                std::fill_n(in_chunk.data() + take, want - take, 0.0f);
            cursor += take;
        }

        stretch.step(want > 0 ? in_chunk.data() : nullptr, pos_pct, out_chunk.data());
        sink(out_chunk.data(), bufsize);
        cursor = clamp_advance(cursor, stretch.skip_after_step(), input.size());
        first = false;
    }
}

std::vector<float> render_channel(
    const std::vector<float> &input,
    const RenderOptions &options,
    const std::vector<Breakpoint> &envelope,
    const ProcessOptions &process_options,
    const std::vector<Breakpoint> &arbitrary_filter) {
    if (input.empty()) return {};

    std::vector<float> output;
    output.reserve(static_cast<std::size_t>(std::ceil(input.size() * options.stretch)) + options.fft_size);

    stream_channel(input, options, envelope, process_options, arbitrary_filter,
                   [&output](const float *data, int frames) {
                       output.insert(output.end(), data, data + frames);
                   });

    return output;
}

// Core streaming loop for stereo. Runs two StreamingStretchers in lockstep and
// delivers each `bufsize` L/R chunk pair to `sink`. Both render_stereo and
// render_stereo_chunked build on this.
void stream_stereo(
    const std::vector<float> &left,
    const std::vector<float> &right,
    const RenderOptions &options,
    const std::vector<Breakpoint> &envelope,
    const ProcessOptions &process_options,
    const std::vector<Breakpoint> &arbitrary_filter,
    const OfflineRenderer::StereoChunkSink &sink) {
    StreamingStretcher stretch_left(options);
    StreamingStretcher stretch_right(options);
    if (!envelope.empty()) {
        stretch_left.set_stretch_envelope(envelope);
        stretch_right.set_stretch_envelope(envelope);
    }
    stretch_left.set_process_options(process_options);
    stretch_right.set_process_options(process_options);
    if (!arbitrary_filter.empty()) {
        stretch_left.set_arbitrary_filter(arbitrary_filter);
        stretch_right.set_arbitrary_filter(arbitrary_filter);
    }

    const int bufsize = stretch_left.bufsize();
    std::vector<float> in_l(stretch_left.max_input_chunk(), 0.0f);
    std::vector<float> in_r(stretch_right.max_input_chunk(), 0.0f);
    std::vector<float> out_l(bufsize, 0.0f);
    std::vector<float> out_r(bufsize, 0.0f);

    bool first = true;
    std::size_t cursor = 0;

    while (true) {
        const float pos_pct = 100.0f * static_cast<float>(cursor) / static_cast<float>(left.size());
        const int want = first ? stretch_left.max_input_chunk() : stretch_left.next_input_size();

        if (want > 0 && cursor >= left.size() && !first) break;

        if (want > 0) {
            const std::size_t avail = left.size() - cursor;
            const std::size_t take = std::min<std::size_t>(avail, static_cast<std::size_t>(want));
            std::copy_n(left.data() + cursor, take, in_l.data());
            std::copy_n(right.data() + cursor, take, in_r.data());
            if (take < static_cast<std::size_t>(want)) {
                std::fill_n(in_l.data() + take, want - take, 0.0f);
                std::fill_n(in_r.data() + take, want - take, 0.0f);
            }
            cursor += take;
        }

        // NB: we drop the per-channel onset coordination the previous offline
        // path did (combining onset_l/onset_r with std::max and feeding both
        // sides) because the StreamingStretcher's step() applies its own onset
        // internally. For independent stereo channels with onset detection on,
        // this means each side reacts to its own transients — usually the more
        // correct behavior anyway.
        stretch_left.step(want > 0 ? in_l.data() : nullptr, pos_pct, out_l.data());
        stretch_right.step(want > 0 ? in_r.data() : nullptr, pos_pct, out_r.data());

        sink(out_l.data(), out_r.data(), bufsize);

        cursor = clamp_advance(cursor, stretch_left.skip_after_step(), left.size());
        first = false;
    }
}

} // anonymous namespace

// ── Public API ─────────────────────────────────────────────────────────────

std::string fft_backend_name() {
#if defined(PAULSTRETCH_USE_PFFFT)
    return "PFFFT";
#elif defined(PAULSTRETCH_USE_ACCELERATE)
    return "Accelerate";
#elif defined(KISSFFT)
    return "KISSFFT";
#else
    return "FFTW";
#endif
}

std::string fft_simd_arch() {
#if defined(PAULSTRETCH_USE_PFFFT)
#if defined(__wasm_simd128__)
    return "WASM_SIMD128";
#else
    return pffft_simd_arch();
#endif
#elif defined(PAULSTRETCH_USE_ACCELERATE)
#if defined(__arm64__) || defined(__aarch64__)
    return "NEON";
#else
    return "SSE";
#endif
#else
    return "scalar";
#endif
}

int fft_simd_size() {
#if defined(PAULSTRETCH_USE_PFFFT)
    return pffft_simd_size();
#elif defined(PAULSTRETCH_USE_ACCELERATE)
    return 4;
#else
    return 1;
#endif
}

namespace {

class AllPass {
public:
    AllPass(float coef = 0.5f) { set(coef); }

    void set(float coef) { a_ = coef * coef; }

    void reset() {
        in1_ = 0.0f;
        in2_ = 0.0f;
        out1_ = 0.0f;
        out2_ = 0.0f;
    }

    float process(float in) {
        const float out = a_ * (in + out2_) - in2_;
        in2_ = in1_;
        in1_ = in;
        out2_ = out1_;
        out1_ = out;
        return out;
    }

private:
    float in1_ = 0.0f;
    float in2_ = 0.0f;
    float out1_ = 0.0f;
    float out2_ = 0.0f;
    float a_ = 0.25f;
};

class Hilbert {
public:
    Hilbert() { reset(); }

    void reset() {
        static constexpr std::array<float, 4> coef_l{
            0.6923877778065f, 0.9360654322959f, 0.9882295226860f, 0.9987488452737f};
        static constexpr std::array<float, 4> coef_r{
            0.4021921162426f, 0.8561710882420f, 0.9722909545651f, 0.9952884791278f};

        old_l_ = 0.0f;
        for (std::size_t i = 0; i < apl_.size(); i++) {
            apl_[i].set(coef_l[i]);
            apr_[i].set(coef_r[i]);
            apl_[i].reset();
            apr_[i].reset();
        }
    }

    void process(float in, float &out1, float &out2) {
        out1 = old_l_;
        out2 = in;
        for (std::size_t i = 0; i < apl_.size(); i++) {
            out1 = apl_[i].process(out1);
            out2 = apr_[i].process(out2);
        }
        old_l_ = in;
    }

private:
    std::array<AllPass, 4> apl_{};
    std::array<AllPass, 4> apr_{};
    float old_l_ = 0.0f;
};

} // anonymous namespace

struct BinauralBeatsProcessor::Impl {
    explicit Impl(float sr) : sample_rate(sr > 0.0f ? sr : 44100.0f) {}

    float beat_frequency(float position_pct) const {
        if (!frequency_envelope.empty())
            return envelope_at(frequency_envelope, position_pct / 100.0f);
        return options.beat_frequency_hz;
    }

    float sample_rate;
    BinauralBeatsOptions options;
    std::vector<Breakpoint> frequency_envelope;
    float hilbert_t = 0.0f;
    Hilbert left_hilbert;
    Hilbert right_hilbert;
};

BinauralBeatsProcessor::BinauralBeatsProcessor(float sample_rate)
    : impl_(std::make_unique<Impl>(sample_rate)) {}

BinauralBeatsProcessor::~BinauralBeatsProcessor() = default;

void BinauralBeatsProcessor::set_options(BinauralBeatsOptions options) {
    options.mono = std::clamp(options.mono, 0.0f, 1.0f);
    options.beat_frequency_hz = std::max(0.0f, options.beat_frequency_hz);
    impl_->options = options;
}

const BinauralBeatsOptions &BinauralBeatsProcessor::options() const {
    return impl_->options;
}

void BinauralBeatsProcessor::set_frequency_envelope(std::vector<Breakpoint> envelope) {
    sort_breakpoints(envelope);
    impl_->frequency_envelope = std::move(envelope);
}

void BinauralBeatsProcessor::clear_frequency_envelope() {
    impl_->frequency_envelope.clear();
}

const std::vector<Breakpoint> &BinauralBeatsProcessor::frequency_envelope() const {
    return impl_->frequency_envelope;
}

void BinauralBeatsProcessor::process(float *left, float *right, int nframes, float position_pct) {
    if (!impl_->options.enabled || !left || !right || nframes <= 0) return;

    const float mono = impl_->options.mono * 0.5f;
    for (int i = 0; i < nframes; i++) {
        const float in_l = left[i];
        const float in_r = right[i];
        left[i] = in_l * (1.0f - mono) + in_r * mono;
        right[i] = in_r * (1.0f - mono) + in_l * mono;
    }

    const float freq = std::max(0.0f, impl_->beat_frequency(position_pct)) * 0.5f;
    for (int i = 0; i < nframes; i++) {
        impl_->hilbert_t = std::fmod(impl_->hilbert_t + freq / impl_->sample_rate, 1.0f);
        const float x = impl_->hilbert_t * 2.0f * kPi;
        const float c = std::cos(x);
        const float s = std::sin(x);

        float h1 = 0.0f;
        float h2 = 0.0f;
        impl_->left_hilbert.process(left[i], h1, h2);
        const float out_l1 = h1 * c + h2 * s;
        const float out_l2 = h1 * c - h2 * s;

        h1 = 0.0f;
        h2 = 0.0f;
        impl_->right_hilbert.process(right[i], h1, h2);
        const float out_r1 = h1 * c - h2 * s;
        const float out_r2 = h1 * c + h2 * s;

        switch (impl_->options.stereo_mode) {
        case BinauralStereoMode::LeftRight:
            left[i] = out_l2;
            right[i] = out_r2;
            break;
        case BinauralStereoMode::RightLeft:
            left[i] = out_l1;
            right[i] = out_r1;
            break;
        case BinauralStereoMode::Symmetric:
            left[i] = (out_l1 + out_r1) * 0.5f;
            right[i] = (out_l2 + out_r2) * 0.5f;
            break;
        }
    }
}

void BinauralBeatsProcessor::reset() {
    impl_->hilbert_t = 0.0f;
    impl_->left_hilbert.reset();
    impl_->right_hilbert.reset();
}

OfflineRenderer::OfflineRenderer(RenderOptions options)
    : options_(sanitize_options(options)) {}

void OfflineRenderer::set_stretch_envelope(std::vector<Breakpoint> envelope) {
    sort_breakpoints(envelope);
    envelope_ = std::move(envelope);
}

void OfflineRenderer::clear_stretch_envelope() {
    envelope_.clear();
}

const std::vector<Breakpoint> &OfflineRenderer::stretch_envelope() const {
    return envelope_;
}

void OfflineRenderer::set_process_options(ProcessOptions options) {
    process_options_ = options;
}

const ProcessOptions &OfflineRenderer::process_options() const {
    return process_options_;
}

void OfflineRenderer::set_arbitrary_filter(std::vector<Breakpoint> filter) {
    sort_breakpoints(filter);
    arbitrary_filter_ = std::move(filter);
}

void OfflineRenderer::clear_arbitrary_filter() {
    arbitrary_filter_.clear();
}

const std::vector<Breakpoint> &OfflineRenderer::arbitrary_filter() const {
    return arbitrary_filter_;
}

const RenderOptions &OfflineRenderer::options() const {
    return options_;
}

std::vector<float> OfflineRenderer::render_mono(const std::vector<float> &input) const {
    return render_channel(input, options_, envelope_, process_options_, arbitrary_filter_);
}

StereoBuffer OfflineRenderer::render_stereo(const std::vector<float> &left, const std::vector<float> &right) const {
    if (left.size() != right.size()) throw std::invalid_argument("left and right channel lengths must match");
    if (left.empty()) return {};

    StereoBuffer output;
    const std::size_t reserve = estimate_output_frames(left.size());
    output.left.reserve(reserve);
    output.right.reserve(reserve);

    stream_stereo(left, right, options_, envelope_, process_options_, arbitrary_filter_,
                  [&output](const float *l, const float *r, int frames) {
                      output.left.insert(output.left.end(), l, l + frames);
                      output.right.insert(output.right.end(), r, r + frames);
                  });

    return output;
}

void OfflineRenderer::render_mono_chunked(const std::vector<float> &input, const ChunkSink &sink) const {
    stream_channel(input, options_, envelope_, process_options_, arbitrary_filter_, sink);
}

void OfflineRenderer::render_stereo_chunked(const std::vector<float> &left,
                                            const std::vector<float> &right,
                                            const StereoChunkSink &sink) const {
    if (left.size() != right.size()) throw std::invalid_argument("left and right channel lengths must match");
    if (left.empty()) return;
    stream_stereo(left, right, options_, envelope_, process_options_, arbitrary_filter_, sink);
}

std::size_t OfflineRenderer::estimate_output_frames(std::size_t input_frames) const {
    return static_cast<std::size_t>(std::ceil(input_frames * options_.stretch))
         + static_cast<std::size_t>(options_.fft_size);
}

} // namespace paulstretch
