/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_RNG_HPP
#define GA_UTILITY_RNG_HPP

#include "utility.hpp"
#include "concepts.hpp"
#include <algorithm>
#include <vector>
#include <span>
#include <random>
#include <limits>
#include <cstdint>
#include <cstddef>
#include <concepts>
#include <atomic>

#ifndef GA_SEED
#define GA_SEED 0x3da99432ab975d26LL
#endif

/** Contains the PRNG classes and functions used for generating random numbers. */
namespace genetic_algorithm::rng
{
    /**
    * Splitmix64 pseudo-random number generator based on https://prng.di.unimi.it/splitmix64.c \n
    * All of the member functions are thread-safe.
    */
    class AtomicSplitmix64
    {
    public:
        using result_type = std::uint64_t;  /**< The generator generates 64 bit integers. */
        using state_type  = std::uint64_t;  /**< The generator has a 64 bit state. */

        /**
        * Create a new generator initialized from a 64 bit seed value.
        * 
        * Instead of creating new instances of this generator, the global prng
        * instance should be used.
        * 
        * @param seed The seed used to initialize the state of the generator.
        */
        explicit constexpr AtomicSplitmix64(state_type seed) noexcept :
            state_(seed)
        {}

        /** Generate the next pseudo-random number of the sequence. Thread-safe. */
        result_type operator()() noexcept;

        /** @returns The smallest value that can be generated. */
        static constexpr result_type min() noexcept;

        /** @returns The largest value that can be generated. */
        static constexpr result_type max() noexcept;

        /** Compare the internal state of 2 generators. @returns True if they are the same. */
        friend constexpr bool operator==(const AtomicSplitmix64&, const AtomicSplitmix64&) = default;

    private:
        std::atomic<state_type> state_;
    };

    /** The pseudo-random number generator class used in the algorithms. */
    using PRNG = AtomicSplitmix64;

    /** The global pseudo-random number generator instance used in the algorithms. */
    inline constinit PRNG prng{ GA_SEED };


    /** Generate a random boolean value from a uniform distribution. */
    inline bool randomBool() noexcept;

    /** Generate a random integer from a uniform distribution on the closed interval [lbound, ubound]. */
    template<std::integral IntType = int>
    inline IntType randomInt(IntType lbound, IntType ubound);

    /** Generate a random floating-point value from a uniform distribution on the closed interval [0.0, 1.0]. */
    template<std::floating_point RealType = double>
    inline RealType randomReal();

    /** Generate a random floating-point value of from a uniform distribution on the closed interval [lbound, ubound]. */
    template<std::floating_point RealType = double>
    inline RealType randomReal(RealType lbound, RealType ubound);

    /** Generate a random floating-point value from a standard normal distribution. */
    template<std::floating_point RealType = double>
    inline RealType randomNormal();

    /** Generate a random floating-point value from a normal distribution with the specified mean and std deviation. */
    template<std::floating_point RealType = double>
    inline RealType randomNormal(RealType mean, RealType SD);

    /** Generate a random integer value from a binomial distribution with the parameters n and p. */
    template<std::integral IntType = int>
    inline IntType randomBinomial(IntType n, double p);

    /** Generate a random index for a container. */
    template<detail::IndexableContainer T>
    inline size_t randomIdx(const T& cont);

    /** Pick a random element from a range. */
    template<std::input_iterator Iter>
    inline Iter randomElement(Iter first, Iter last);

    /** Generate @p count unique integers from the half-open range [@p lbound, @p ubound). */
    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count);

    /** Select an index based on a discrete CDF. */
    template<std::floating_point T>
    inline size_t sampleCdf(std::span<const T> cdf);

} // namespace genetic_algorithm::rng


/* IMPLEMENTATION */

#include <utility>
#include <numeric>
#include <iterator>
#include <cmath>
#include <climits>

namespace genetic_algorithm::rng
{
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
        constexpr size_t nbits = CHAR_BIT * sizeof(PRNG::result_type);
        constexpr auto msb_mask = PRNG::result_type{ 1 } << (nbits - 1);

        return static_cast<bool>(rng::prng() & msb_mask);
    }

    template<std::integral IntType>
    IntType randomInt(IntType lbound, IntType ubound)
    {
        GA_ASSERT(lbound <= ubound);

        return std::uniform_int_distribution{ lbound, ubound }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal()
    {
        return std::uniform_real_distribution<RealType>{ 0.0, std::nextafter(1.0, 2.0) }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal(RealType lbound, RealType ubound)
    {
        GA_ASSERT(lbound <= ubound);

        ubound = std::nextafter(ubound, std::numeric_limits<RealType>::max());

        return std::uniform_real_distribution{ lbound, ubound }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal()
    {
        thread_local std::normal_distribution<RealType> dist{ 0.0, 1.0 }; // keep the distribution for the state

        return dist(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal(RealType mean, RealType SD)
    {
        GA_ASSERT(SD >= 0.0);

        return (SD == 0.0) ? mean : std::normal_distribution{ mean, SD }(rng::prng);
    }

    template<std::integral IntType>
    IntType randomBinomialApprox(IntType n, double p)
    {
        GA_ASSERT(n >= 0);
        GA_ASSERT(0.0 <= p && p <= 1.0);

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
        GA_ASSERT(n >= 0);
        GA_ASSERT(0.0 <= p && p <= 1.0);

        return std::binomial_distribution{ n, p }(rng::prng);
    }

    template<std::integral IntType>
    IntType randomBinomial(IntType n, double p)
    {
        GA_ASSERT(n >= 0);
        GA_ASSERT(0.0 <= p && p <= 1.0);

        const double mean = n * p;

        return (mean >= 2.0) ?
            rng::randomBinomialApprox(n, p) :
            rng::randomBinomialExact(n, p);
    }

    template<detail::IndexableContainer T>
    size_t randomIdx(const T& container)
    {
        GA_ASSERT(!container.empty());

        return std::uniform_int_distribution{ 0_sz, container.size() - 1 }(rng::prng);
    }

    template<std::input_iterator Iter>
    Iter randomElement(Iter first, Iter last)
    {
        GA_ASSERT(std::distance(first, last) > 0);

        const auto max_offset = std::distance(first, last) - 1;
        const auto offset = rng::randomInt(0_pd, max_offset);

        return std::next(first, offset);
    }

    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count)
    {
        const auto range_len = size_t(ubound - lbound);

        GA_ASSERT(ubound >= lbound);
        GA_ASSERT(range_len >= count);

        std::vector<bool> is_selected(range_len);
        std::vector<IntType> numbers(count);

        for (IntType i = ubound - IntType(count); i < ubound; i++)
        {
            const auto num = rng::randomInt(lbound, i);
            const auto spos = size_t(num - lbound);
            const auto npos = size_t(i + IntType(count) - ubound);

            numbers[npos] = is_selected[spos] ? i : num;
            is_selected[numbers[npos] - lbound] = true;
        }

        return numbers;
    }

    template<std::floating_point T>
    size_t sampleCdf(std::span<const T> cdf)
    {
        GA_ASSERT(!cdf.empty());

        const auto limit = rng::randomReal<T>(0.0, cdf.back()); // use cdf.back() in case it's not exactly 1.0
        const auto selected = std::lower_bound(cdf.begin(), cdf.end(), limit);

        return size_t(selected - cdf.begin());
    }

} // namespace genetic_algorithm::rng

#endif // !GA_UTILITY_RNG_HPP