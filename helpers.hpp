#pragma once
#include <type_traits>
#include <limits>
#include <cstdint>

namespace fp {

// Map any bit count to {8,16,32}
template<int B> struct BucketBits {
    static constexpr int value = (B <= 8) ? 8 : (B <= 16) ? 16 : 32;
};

template<int Bucket> struct StorageFromBucket;
template<> struct StorageFromBucket<8>  { using type = int8_t;  };
template<> struct StorageFromBucket<16> { using type = int16_t; };
template<> struct StorageFromBucket<32> { using type = int32_t; };

template<int B> struct StorageForBits {
    using type = typename StorageFromBucket< BucketBits<B>::value >::type;
};

// Convenient alias for StorageForBits<B>::type
template<int B>
using Storage_t = typename StorageForBits<B>::type;

template<int B, int Want>
struct IsBucket : std::integral_constant<bool, (BucketBits<B>::value == Want)> {};

// Convenient alias for enable_if
template<bool B>
using EnableIf = typename std::enable_if<B, int>::type;

// Saturating cast
template<typename To, typename From>
constexpr To sat_cast(From v) {
    using Lim = std::numeric_limits<To>;
    long long w = static_cast<long long>(v);
    if (w > static_cast<long long>(Lim::max())) return Lim::max();
    if (w < static_cast<long long>(Lim::min())) return Lim::min();
    return static_cast<To>(w);
}

// Priority tag ladder
template<int N> struct priority_tag : priority_tag<N-1> {};
template<> struct priority_tag<0> {};

// signed round-to-nearest, ties away from zero (simple variant)
inline auto round_shift = [](long long x, int s) -> long long {
    if (s <= 0) return (s==0 ? x : (x << (-s)));
    long long bias = 1ll << (s - 1);
    return (x >= 0) ? (x + bias) >> s : (x - bias) >> s;
};

} // namespace fp
