#pragma once
#include <cstdint>
#include <cmath>
#include <cstdio>   // for FILE*, fopen, fprintf, fclose
#include <array>

#ifndef OBSERVER_ENABLE
#define OBSERVER_ENABLE 1
#endif

// NEW: enable histogram support (fixed bins, low RAM)
#ifndef OBSERVER_HIST_ENABLE
#define OBSERVER_HIST_ENABLE 1
#endif

// NEW: histogram config
#ifndef OBSERVER_HIST_BINS
#define OBSERVER_HIST_BINS 64
#endif

#ifndef OBSERVER_HIST_ABSMAX_CLIP   // clip values with |x| > this into last bin
#define OBSERVER_HIST_ABSMAX_CLIP  (1e6)
#endif

namespace obs {

// ----------------- RunningStats -----------------
struct RunningStats {
    uint64_t n = 0;
    double mean = 0.0;
    double m2 = 0.0;
    double minv =  INFINITY;
    double maxv = -INFINITY;

    inline void reset() {
        n = 0; mean = 0.0; m2 = 0.0; minv = INFINITY; maxv = -INFINITY;
    }
    inline void update(double x) {
        if (x < minv) minv = x;
        if (x > maxv) maxv = x;
        n++;
        double d = x - mean;
        mean += d / double(n);
        double d2 = x - mean;
        m2 += d * d2;
    }
    inline double variance() const { return (n > 1) ? (m2 / double(n - 1)) : 0.0; }
    inline double stddev()   const { return std::sqrt(variance()); }
    inline double min()      const { return minv; }
    inline double max()      const { return maxv; }
    inline uint64_t count()  const { return n; }
    inline double absmax()   const { return std::max(std::abs(minv), std::abs(maxv)); }
};

// ----------------- Histogram (optional) -----------------
#if OBSERVER_HIST_ENABLE
struct Histogram {
    // Symmetric linear bins over [-range, +range].
    // Values outside are clipped into the first/last bin.
    double range = 1.0; // half-span (i.e., covers [-range, +range])
    std::array<uint64_t, OBSERVER_HIST_BINS> bins{};

    inline void reset(double new_range = 1.0) {
        range = (new_range > 1e-30) ? new_range : 1.0;
        for (auto &b : bins) b = 0;
    }
    inline void update(double x) {
        double r = range;
        if (r <= 0) r = 1.0;
        // clip
        if (x <= -OBSERVER_HIST_ABSMAX_CLIP) x = -OBSERVER_HIST_ABSMAX_CLIP;
        if (x >=  OBSERVER_HIST_ABSMAX_CLIP) x =  OBSERVER_HIST_ABSMAX_CLIP;

        // map x from [-r, r] -> [0, OBSERVER_HIST_BINS-1]
        double t = (x + r) / (2.0 * r);                 // [0,1]
        long idx = (long)std::floor(t * OBSERVER_HIST_BINS);
        if (idx < 0) idx = 0;
        if (idx >= (long)OBSERVER_HIST_BINS) idx = OBSERVER_HIST_BINS - 1;
        bins[(size_t)idx]++;
    }
    // choose range from absmax with a little headroom
    inline void fit_range_from_absmax(double absmax, double headroom = 1.05) {
        double r = absmax * headroom;
        if (r < 1e-6) r = 1.0; // avoid degenerate
        range = r;
    }
};
#endif

// ----------------- Ops & Channel -----------------
enum class Op : uint8_t { Value=0, Sum, Diff, Prod, Quot, Accum, _COUNT };

template<typename Tag>
struct Channel {
#if OBSERVER_ENABLE
    struct OpBucket {
        RunningStats stats;
    #if OBSERVER_HIST_ENABLE
        Histogram hist;
    #endif
        inline void reset() {
            stats.reset();
        #if OBSERVER_HIST_ENABLE
            hist.reset(hist.range);
        #endif
        }
        inline void update(double x) {
            stats.update(x);
        #if OBSERVER_HIST_ENABLE
            hist.update(x);
        #endif
        }
    };

    static inline OpBucket buckets[(size_t)Op::_COUNT]{};
    static inline const char*  name = Tag::name;

    static inline void log(Op op, double v) {
        buckets[(size_t)op].update(v);
    }
    static inline const RunningStats& get(Op op) {
        return buckets[(size_t)op].stats;
    }
    static inline void reset_all() {
        for (size_t i=0; i<(size_t)Op::_COUNT; ++i) buckets[i].reset();
    }

    // Prepare hist ranges after a warmup pass (optional).
    // Call this if you want histograms scaled to observed ranges before a second pass.
    static inline void fit_hist_ranges(double headroom = 1.05) {
    #if OBSERVER_HIST_ENABLE
        for (size_t i=0; i<(size_t)Op::_COUNT; ++i) {
            double amax = buckets[i].stats.absmax();
            buckets[i].hist.fit_range_from_absmax(amax, headroom);
        }
    #else
        (void)headroom;
    #endif
    }

    // ------------- CSV dump -------------
    // Writes one CSV per op to `basepath + "." + opname + ".csv"`
    // Format:
    //   header: key,value
    //   stats rows (min,max,mean,stddev,count,absmax,range)
    //   blank line
    //   header: bin_left,bin_right,count
    //   per-bin rows
    static inline bool dump_csv(const char* basepath) {
        static const char* kOpName[] = {"value","sum","diff","prod","quot","accum"};
        bool ok_all = true;

        for (size_t i=0; i<(size_t)Op::_COUNT; ++i) {
            char path[256];
            std::snprintf(path, sizeof(path), "%s.%s.csv", basepath, kOpName[i]);
            FILE* f = std::fopen(path, "wb");
            if (!f) { ok_all = false; continue; }

            const auto& rs = buckets[i].stats;

            // stats section
            std::fprintf(f, "key,value\n");
            std::fprintf(f, "tag,%s\n", name ? name : "unknown");
            std::fprintf(f, "op,%s\n", kOpName[i]);
            std::fprintf(f, "count,%llu\n", (unsigned long long)rs.count());
            std::fprintf(f, "min,%.17g\n", rs.min());
            std::fprintf(f, "max,%.17g\n", rs.max());
            std::fprintf(f, "mean,%.17g\n", rs.mean);
            std::fprintf(f, "stddev,%.17g\n", rs.stddev());
            std::fprintf(f, "absmax,%.17g\n", rs.absmax());

        #if OBSERVER_HIST_ENABLE
            const auto& h = buckets[i].hist;
            std::fprintf(f, "hist_range,%.17g\n", h.range);
            std::fprintf(f, "\n"); // spacer

            // histogram section
            std::fprintf(f, "bin_left,bin_right,count\n");
            const double r  = h.range;
            const double bw = (2.0 * r) / double(OBSERVER_HIST_BINS);
            for (int b = 0; b < OBSERVER_HIST_BINS; ++b) {
                double left  = -r + b * bw;
                double right = left + bw;
                std::fprintf(f, "%.17g,%.17g,%llu\n",
                             left, right,
                             (unsigned long long)h.bins[(size_t)b]);
            }
        #else
            std::fprintf(f, "\n");
            std::fprintf(f, "bin_left,bin_right,count\n"); // empty hist
        #endif
            std::fclose(f);
        }
        return ok_all;
    }

#else
    static inline void log(Op, double) {}
    static inline bool dump_csv(const char*) { return true; }
    static inline void reset_all() {}
    static inline void fit_hist_ranges(double = 1.05) {}
    static inline const RunningStats& get(Op) {
        static RunningStats dummy; return dummy;
    }
#endif
};

// ---- QSuggest unchanged from previous message ----
struct QSuggest {
    int total_bits;
    bool signedness;
    inline std::pair<int,int> suggest_from_absmax(double absmax) const {
        if (absmax < 1e-30) absmax = 0.0;
        int sign_bits = signedness ? 1 : 0;
        double eps = 1e-30;
        int int_bits = (absmax <= 1.0) ? 1 : int(std::ceil(std::log2(absmax + eps))) + 1;
        int frac_bits = total_bits - sign_bits - int_bits;
        if (frac_bits < 0) { frac_bits = 0; int_bits = total_bits - sign_bits; }
        return {int_bits, frac_bits};
    }
};

// ---- Helper functions for printing statistics ----

// Get operation name as string
inline const char* op_name(Op op) {
    static const char* names[] = {"Value", "Sum", "Diff", "Prod", "Quot", "Accum"};
    return names[(size_t)op];
}

// Print statistics for a single operation on a channel
template<typename Tag>
inline void print_channel_op(Op op, const QSuggest& qs16, const QSuggest& qs32, bool add_newline = false) {
#if OBSERVER_ENABLE
    const auto& stats = Channel<Tag>::get(op);
    auto [i16, f16] = qs16.suggest_from_absmax(stats.absmax());
    auto [i32, f32] = qs32.suggest_from_absmax(stats.absmax());

    std::printf("[%s/%s] count=%llu min=%.6g max=%.6g mean=%.6g std=%.6g | Q16: Q%d.%d | Q32: Q%d.%d\n%s",
               Tag::name, op_name(op),
               (unsigned long long)stats.count(), stats.min(), stats.max(), stats.mean, stats.stddev(),
               i16, f16, i32, f32,
               add_newline ? "\n" : "");
#else
    (void)op; (void)qs16; (void)qs32; (void)add_newline;
#endif
}

// Print statistics for multiple operations on a channel
template<typename Tag>
inline void print_channel_ops(std::initializer_list<Op> ops, const QSuggest& qs16, const QSuggest& qs32) {
#if OBSERVER_ENABLE
    auto it = ops.begin();
    auto end = ops.end();
    while (it != end) {
        auto next_it = it;
        ++next_it;
        bool is_last = (next_it == end);
        print_channel_op<Tag>(*it, qs16, qs32, is_last);
        it = next_it;
    }
#else
    (void)ops; (void)qs16; (void)qs32;
#endif
}

} // namespace obs

// -------- Macros as in previous message (unchanged) --------
#if OBSERVER_ENABLE
    #define OBS_LOG(TAG, OP, VALUE) ::obs::Channel<TAG>::log(::obs::Op::OP, double(VALUE))
#else
    #define OBS_LOG(TAG, OP, VALUE) do{}while(0)
#endif

#define OBS_VAL(TAG, x)         ([&](){ auto _v=(x); OBS_LOG(TAG, Value, _v); return _v; }())
#define OBS_ADD(TAG, x,y)       ([&](){ auto _a=(x), _b=(y); auto _r=_a+_b; OBS_LOG(TAG, Sum, _r); return _r; }())
#define OBS_SUB(TAG, x,y)       ([&](){ auto _a=(x), _b=(y); auto _r=_a-_b; OBS_LOG(TAG, Diff,_r); return _r; }())
#define OBS_MUL(TAG, x,y)       ([&](){ auto _a=(x), _b=(y); auto _r=_a*_b; OBS_LOG(TAG, Prod,_r); return _r; }())
#define OBS_DIV(TAG, x,y)       ([&](){ auto _a=(x), _b=(y); auto _r=_a/_b; OBS_LOG(TAG, Quot,_r); return _r; }())
#define OBS_ACCUM(TAG, acc,term)([&](){ auto _t=(term); (acc)+=_t; OBS_LOG(TAG, Accum,(acc)); return (acc); }())

namespace obs {
template<typename T, typename Tag>
struct Observed {
    T v;
    explicit Observed(T x) : v(x) { Channel<Tag>::log(Op::Value, double(x)); }
    operator T() const { return v; }
    friend Observed operator+(Observed a, Observed b){ T r=a.v+b.v; Channel<Tag>::log(Op::Sum,  r); return Observed(r); }
    friend Observed operator-(Observed a, Observed b){ T r=a.v-b.v; Channel<Tag>::log(Op::Diff, r); return Observed(r); }
    friend Observed operator*(Observed a, Observed b){ T r=a.v*b.v; Channel<Tag>::log(Op::Prod, r); return Observed(r); }
    friend Observed operator/(Observed a, Observed b){ T r=a.v/b.v; Channel<Tag>::log(Op::Quot, r); return Observed(r); }
    Observed& operator+=(Observed rhs){ v += rhs.v; Channel<Tag>::log(Op::Accum, v); return *this; }
    Observed& operator-=(Observed rhs){ v -= rhs.v; Channel<Tag>::log(Op::Diff,  v); return *this; }
    Observed& operator*=(Observed rhs){ v *= rhs.v; Channel<Tag>::log(Op::Prod,  v); return *this; }
    Observed& operator/=(Observed rhs){ v /= rhs.v; Channel<Tag>::log(Op::Quot,  v); return *this; }
};
} // namespace obs

// Example Tag definitions
/*
 struct DotTag { static constexpr const char* name = "dot"; };

float AfcNlmsNode::dot(float* a, float* b, uint32_t n) {
    float acc = 0.0f;
    for (uint32_t i=0; i<n; ++i) {
        float prod = OBS_MUL(DotTag, a[i], b[i]);
        OBS_ACCUM(DotTag, acc, prod);
    }
    return acc;
}

// After a warmup pass (optional), set per-op histogram ranges from observed absmax:
obs::Channel<DotTag>::fit_hist_ranges(1.10); // 10% headroom

// After your real run(s), dump CSVs like: dot.value.csv, dot.prod.csv, dot.accum.csv, ...
obs::Channel<DotTag>::dump_csv("dot");

// Suggest Q from stats:
const auto& accum = obs::Channel<DotTag>::get(obs::Op::Accum);
obs::QSuggest qs{16, true};
auto [i_bits, f_bits] = qs.suggest_from_absmax(accum.absmax());

 */