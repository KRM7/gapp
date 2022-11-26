/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_RNG_HPP
#define GA_UTILITY_RNG_HPP

#include "utility.hpp"
#include "concepts.hpp"
#include <vector>
#include <random>
#include <limits>
#include <cstdint>
#include <cstddef>
#include <concepts>
#include <atomic>

#ifndef GA_SEED
#define GA_SEED 0x3da99432ab975d26LL
#endif

/** Contains the PRNG classes and functions for generating random numbers. */
namespace genetic_algorithm::rng
{
    /**  Splitmix64 PRNG adapted from https://prng.di.unimi.it/splitmix64.c */
    class AtomicSplitmix64
    {
    public:
        using result_type = std::uint64_t;
        using state_type  = std::uint64_t;

        explicit constexpr AtomicSplitmix64(state_type seed) noexcept;

        result_type operator()() noexcept;

        static constexpr result_type min() noexcept;
        static constexpr result_type max() noexcept;

    private:
        std::atomic<state_type> state_;
    };

    /** The PRNG type used in the genetic algorithms. */
    using PRNG = AtomicSplitmix64;

    /** PRNG instance used in the genetic algorithms. */
    extern PRNG prng;

    /** Generates a random boolean value from a uniform distribution. */
    inline bool randomBool() noexcept;

    /** Generates a random integer of type IntType from a uniform distribution on the closed interval [lbound, ubound]. */
    template<std::integral IntType = int>
    inline IntType randomInt(IntType lbound, IntType ubound);

    /** Generates a random floating-point value of type RealType from a uniform distribution on the interval [0.0, 1.0). */
    template<std::floating_point RealType = double>
    inline RealType randomReal();

    /** Generates a random floating-point value of type RealType from a uniform distribution on the interval [lbound, ubound). */
    template<std::floating_point RealType = double>
    inline RealType randomReal(RealType lbound, RealType ubound);

    /** Generates a random floating-point value of type RealType from a standard normal distribution. */
    template<std::floating_point RealType = double>
    inline RealType randomNormal();

    /** Generates a random floating-point value of type RealType from a normal distribution with the parameters mean and SD. */
    template<std::floating_point RealType = double>
    inline RealType randomNormal(RealType mean, RealType SD);

    /** Generates a random integer from a binomial distribution with the parameters n and p. */
    template<std::integral IntType = int>
    inline IntType randomBinomial(IntType n, double p);

    /** Generates a random idx for a container. */
    template<detail::IndexableContainer T>
    inline size_t randomIdx(const T& cont);

    /** Pick a random element from a range. */
    template<std::input_iterator Iter>
    inline auto randomElement(Iter first, Iter last);

    /** Generates @p count unique integers from the range [lbound, ubound). */
    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count);

    /** Select an index based on a discrete CDF. */
    template<std::floating_point T>
    inline size_t sampleCdf(const std::vector<T>& cdf);

} // namespace genetic_algorithm::rng


/* IMPLEMENTATION */

#include <utility>
#include <numeric>
#include <iterator>
#include <cmath>
#include <climits>
#include <cassert>

namespace genetic_algorithm::rng
{
    constexpr inline AtomicSplitmix64::AtomicSplitmix64(state_type seed) noexcept
        : state_(seed)
    {}

    inline AtomicSplitmix64::result_type AtomicSplitmix64::operator()() noexcept
    {
        result_type z = state_.fetch_add(0x9e3779b97f4a7c15, std::memory_order_relaxed) + 0x9e3779b97f4a7c15;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;

        return z ^ (z >> 31);
    }

    inline constexpr AtomicSplitmix64::result_type AtomicSplitmix64::min() noexcept
    {
        return std::numeric_limits<result_type>::min();
    }

    inline constexpr AtomicSplitmix64::result_type AtomicSplitmix64::max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

    bool randomBool() noexcept
    {
        constexpr size_t nbytes = sizeof(PRNG::result_type);
        constexpr size_t nbits = CHAR_BIT * nbytes;
        constexpr auto msb_mask = PRNG::result_type{ 1 } << (nbits - 1);

        return rng::prng() & msb_mask;
    }

    template<std::integral IntType>
    IntType randomInt(IntType lbound, IntType ubound)
    {
        assert(lbound <= ubound);

        return std::uniform_int_distribution{ lbound, ubound }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal()
    {
        return std::uniform_real_distribution<RealType>{ 0.0, 1.0 }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal(RealType lbound, RealType ubound)
    {
        assert(lbound <= ubound);

        return std::uniform_real_distribution{ lbound, ubound }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal()
    {
        return std::normal_distribution<RealType>{ 0.0, 1.0 }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal(RealType mean, RealType SD)
    {
        assert(SD > 0.0);

        return std::normal_distribution{ mean, SD }(rng::prng);
    }

    template<std::integral IntType>
    IntType randomBinomialApprox(IntType n, double p)
    {
        assert(0.0 <= p && p <= 1.0);

        const double mean = n * p;
        const double SD = std::sqrt(mean * (1.0 - p));

        const double accept_min = -0.5;
        const double accept_max = n + 0.5;

        double rand = rng::randomNormal(mean, SD);
        while (!(accept_min < rand && rand < accept_max))
        {
            rand = rng::randomNormal(mean, SD);
        }

        return static_cast<IntType>(std::round(rand));
    }

    template<std::integral IntType>
    IntType randomBinomialExact(IntType n, double p)
    {
        return std::binomial_distribution{ n, p }(rng::prng);
    }

    template<std::integral IntType>
    IntType randomBinomial(IntType n, double p)
    {
        assert(0.0 <= p && p <= 1.0);

        const double mean = n * p;

        return (mean >= 2.0) ?
            rng::randomBinomialApprox(n, p) :
            rng::randomBinomialExact(n, p);
    }

    template<detail::IndexableContainer T>
    size_t randomIdx(const T& container)
    {
        assert(!container.empty());

        return std::uniform_int_distribution{ 0_sz, container.size() - 1 }(rng::prng);
    }

    template<std::input_iterator Iter>
    auto randomElement(Iter first, Iter last)
    {
        assert(std::distance(first, last) > 0);

        const auto max_offset = std::distance(first, last) - 1;
        const auto offset = rng::randomInt(0_pd, max_offset);

        return *std::next(first, offset);
    }

    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count)
    {
        const auto range_len = size_t(ubound - lbound);

        assert(ubound >= lbound);
        assert(range_len >= count);

        std::vector is_selected(range_len, false);
        std::vector numbers(count, ubound);

        for (IntType i = ubound - IntType(count); i < ubound; i++)
        {
            const auto num = rng::randomInt(lbound, i);
            const auto spos = size_t(num - lbound);
            const auto npos = size_t(i + IntType(count) - ubound);

            numbers[npos] = is_selected[spos] ? i : num;
            is_selected[spos] = true;
        }

        return numbers;
    }

    template<std::floating_point T>
    size_t sampleCdf(const std::vector<T>& cdf)
    {
        assert(!cdf.empty());

        const auto limit = rng::randomReal<T>(0.0, cdf.back()); // use cdf.back() in case it's not exactly 1.0
        const auto selected = std::lower_bound(cdf.begin(), cdf.end(), limit);

        return size_t(selected - cdf.begin());
    }

} // namespace genetic_algorithm::rng

#endif // !GA_UTILITY_RNG_HPP