#include "paulstretch/paulstretch.h"

#include <algorithm>
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
    std::unique_ptr<Stretcher> stretch;
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
    Stretcher &s = *impl_->stretch;
    // Always record position so the envelope is evaluated at the caller's
    // current cursor — important for seek to land on the right envelope
    // value on the very next step.
    const int natural_n = s.get_nsamples(position_pct);
    const int n = impl_->first_step ? s.max_bufsize() : natural_n;
    impl_->first_step = false;

    const float onset = s.process(n > 0 ? input : nullptr, n);
    s.here_is_onset(onset);
    std::copy_n(s.out_buf, s.bufsize(), output);
    impl_->skip_after = s.get_skip_nsamples();
    return onset;
}

void StreamingStretcher::set_stretch_envelope(std::vector<Breakpoint> envelope) {
    std::sort(envelope.begin(), envelope.end(),
              [](const Breakpoint &a, const Breakpoint &b) { return a.position < b.position; });
    impl_->envelope = std::move(envelope);
    // Stretcher holds a raw pointer to the envelope vector. Rebuild so the
    // new envelope is wired in cleanly. (Rebuilding also resets DSP state,
    // which matches the offline behaviour where the renderer constructs a
    // fresh Stretcher per render.)
    impl_->rebuild();
}

void StreamingStretcher::clear_stretch_envelope() {
    impl_->envelope.clear();
    impl_->rebuild();
}

const std::vector<Breakpoint> &StreamingStretcher::stretch_envelope() const {
    return impl_->envelope;
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

std::vector<float> render_channel(
    const std::vector<float> &input,
    const RenderOptions &options,
    const std::vector<Breakpoint> &envelope) {
    if (input.empty()) return {};

    StreamingStretcher stretch(options);
    if (!envelope.empty()) stretch.set_stretch_envelope(envelope);

    const int bufsize = stretch.bufsize();
    std::vector<float> out_chunk(bufsize, 0.0f);
    std::vector<float> in_chunk(stretch.max_input_chunk(), 0.0f);
    std::vector<float> output;
    output.reserve(static_cast<std::size_t>(std::ceil(input.size() * options.stretch)) + bufsize);

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
        output.insert(output.end(), out_chunk.begin(), out_chunk.end());
        cursor = clamp_advance(cursor, stretch.skip_after_step(), input.size());
        first = false;
    }

    return output;
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

OfflineRenderer::OfflineRenderer(RenderOptions options)
    : options_(sanitize_options(options)) {}

void OfflineRenderer::set_stretch_envelope(std::vector<Breakpoint> envelope) {
    std::sort(envelope.begin(), envelope.end(),
              [](const Breakpoint &a, const Breakpoint &b) { return a.position < b.position; });
    envelope_ = std::move(envelope);
}

void OfflineRenderer::clear_stretch_envelope() {
    envelope_.clear();
}

const std::vector<Breakpoint> &OfflineRenderer::stretch_envelope() const {
    return envelope_;
}

const RenderOptions &OfflineRenderer::options() const {
    return options_;
}

std::vector<float> OfflineRenderer::render_mono(const std::vector<float> &input) const {
    return render_channel(input, options_, envelope_);
}

StereoBuffer OfflineRenderer::render_stereo(const std::vector<float> &left, const std::vector<float> &right) const {
    if (left.size() != right.size()) throw std::invalid_argument("left and right channel lengths must match");
    if (left.empty()) return {};

    StreamingStretcher stretch_left(options_);
    StreamingStretcher stretch_right(options_);
    if (!envelope_.empty()) {
        stretch_left.set_stretch_envelope(envelope_);
        stretch_right.set_stretch_envelope(envelope_);
    }

    const int bufsize = stretch_left.bufsize();
    std::vector<float> in_l(stretch_left.max_input_chunk(), 0.0f);
    std::vector<float> in_r(stretch_right.max_input_chunk(), 0.0f);
    std::vector<float> out_l(bufsize, 0.0f);
    std::vector<float> out_r(bufsize, 0.0f);

    StereoBuffer output;
    const std::size_t reserve = estimate_output_frames(left.size());
    output.left.reserve(reserve);
    output.right.reserve(reserve);

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

        output.left.insert(output.left.end(), out_l.begin(), out_l.end());
        output.right.insert(output.right.end(), out_r.begin(), out_r.end());

        cursor = clamp_advance(cursor, stretch_left.skip_after_step(), left.size());
        first = false;
    }

    return output;
}

std::size_t OfflineRenderer::estimate_output_frames(std::size_t input_frames) const {
    return static_cast<std::size_t>(std::ceil(input_frames * options_.stretch))
         + static_cast<std::size_t>(options_.fft_size);
}

} // namespace paulstretch
