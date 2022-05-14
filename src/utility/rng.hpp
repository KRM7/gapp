/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_RANDOM_HPP
#define GA_RANDOM_HPP

#include "utility.hpp"
#include "concepts.hpp"
#include <vector>
#include <random>
#include <iterator>
#include <cstdint>
#include <cstddef>
#include <concepts>
#include <atomic>

/** Contains the PRNG classes and functions for generating random numbers. */
namespace genetic_algorithm::rng
{
    /**  Splitmix64 PRNG adapted from https://prng.di.unimi.it/splitmix64.c */
    class Splitmix64
    {
    public:
        using result_type = uint64_t;
        using state_type = uint64_t;

        explicit constexpr Splitmix64(state_type seed) noexcept;

        result_type operator()() noexcept;

        static constexpr result_type min() noexcept;
        static constexpr result_type max() noexcept;

    private:
        std::atomic<state_type> state_;
    };

    /** The PRNG type used in the genetic algorithms. */
    using PRNG = Splitmix64;

    /** PRNG instance used in the genetic algorithms. */
    inline PRNG prng{ std::random_device{}() };

    /** Generates a random floating-point value of type RealType from a uniform distribution on the interval [0.0, 1.0). */
    template<std::floating_point RealType = double>
    RealType randomReal();

    /** Generates a random floating-point value of type RealType from a uniform distribution on the interval [l_bound, u_bound). */
    template<std::floating_point RealType = double>
    RealType randomReal(RealType l_bound, RealType u_bound);

    /** Generates a random floating-point value of type RealType from a standard normal distribution. */
    template<std::floating_point RealType = double>
    RealType randomNormal();

    /** Generates a random floating-point value of type RealType from a normal distribution with the parameters mean and SD. */
    template<std::floating_point RealType = double>
    RealType randomNormal(RealType mean, RealType SD);

    /** Generates a random integer from a binomial distribution with the parameters n and p. */
    template<std::integral IntType = int>
    IntType randomBinomial(IntType n, double p);

    /** Generates a random integer from a distribution that is the approximation of a binomial distribution with the parameters n and p. */
    template<std::integral IntType = int>
    IntType randomBinomialApprox(IntType n, double p);

    /** Generates a random integer of type IntType from a uniform distribution on the closed interval [l_bound, u_bound]. */
    template<std::integral IntType = int>
    IntType randomInt(IntType l_bound, IntType u_bound);

    /** Generates a random idx for a container. */
    template<detail::IndexableContainer T>
    size_t randomIdx(const T& cont);

    /** Pick a random element from a container. */
    template<detail::Container T>
    auto randomElement(const T& cont) -> typename T::value_type;

    /** Pick a random element from a range. */
    template<std::input_iterator Iter>
    auto randomElement(Iter first, Iter last);

    /** Generates a random boolean value from a uniform distribution. */
    inline bool randomBool() noexcept;

    /** Generates n unique integers from the range [l_bound, u_bound). */
    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType l_bound, IntType u_bound, size_t n);

    /** Sample a point from a uniform distribution on a unit simplex in dim dimensions. */
    inline std::vector<double> randomSimplexPoint(size_t dim);

    /** Select an index based on a discrete CDF. */
    inline size_t sampleCdf(const std::vector<double>& cdf);

} // namespace genetic_algorithm::rng


/* IMPLEMENTATION */

#include <limits>
#include <utility>
#include <numeric>
#include <mutex>
#include <climits>
#include <cassert>

namespace genetic_algorithm::rng
{
    constexpr inline Splitmix64::Splitmix64(state_type seed) noexcept
        : state_(seed)
    {
    }

    inline Splitmix64::result_type Splitmix64::operator()() noexcept
    {
        result_type z = state_.fetch_add(0x9e3779b97f4a7c15, std::memory_order::acquire) + 0x9e3779b97f4a7c15;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;

        return z ^ (z >> 31);
    }

    inline constexpr Splitmix64::result_type Splitmix64::min() noexcept
    {
        return std::numeric_limits<result_type>::min();
    }

    inline constexpr Splitmix64::result_type Splitmix64::max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

    template<std::floating_point RealType>
    RealType randomReal()
    {
        return std::uniform_real_distribution<RealType>{ 0.0, 1.0 }(prng);
    }

    template<std::floating_point RealType>
    RealType randomReal(RealType l_bound, RealType u_bound)
    {
        assert(l_bound <= u_bound);

        return std::uniform_real_distribution{ l_bound, u_bound }(prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal()
    {
        return std::normal_distribution<RealType>{ 0.0, 1.0 }(prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal(RealType mean, RealType SD)
    {
        assert(SD > 0.0);

        return std::normal_distribution{ mean, SD }(prng);
    }

    template<std::integral IntType>
    IntType randomBinomial(IntType n, double p)
    {
        return std::binomial_distribution{ n, p }(prng);
    }

    template<std::integral IntType>
    IntType randomBinomialApprox(IntType n, double p)
    {
        double mean = n * p;
        if (mean >= 2.0)
        {
            double SD = std::sqrt(mean * (1.0 - p));

            double r = randomNormal(mean, SD);
            while (r <= -0.5) { r = randomNormal(mean, SD); }

            IntType result = IntType(std::round(r));
            return std::min(result, n);
        }
        else
        {
            return randomBinomial(n, p);
        }
    }

    template<std::integral IntType>
    IntType randomInt(IntType l_bound, IntType u_bound)
    {
        assert(l_bound <= u_bound);

        return std::uniform_int_distribution{ l_bound, u_bound }(prng);
    }

    template<detail::IndexableContainer T>
    size_t randomIdx(const T& container)
    {
        assert(!container.empty());

        return std::uniform_int_distribution{ 0_sz, container.size() - 1 }(prng);
    }

    template<detail::Container T>
    auto randomElement(const T& cont) -> typename T::value_type
    {
        size_t n = rng::randomInt(0_sz, cont.size() - 1);

        return *std::next(cont.begin(), n);
    }

    template<std::input_iterator Iter>
    auto randomElement(Iter first, Iter last)
    {
        auto n = rng::randomInt<ptrdiff_t>(0, std::distance(first, last) - 1);

        return *std::next(first, n);
    }

    bool randomBool() noexcept
    {
        return prng() & 1;
    }

    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType l_bound, IntType u_bound, size_t n)
    {
        assert(l_bound <= u_bound);
        assert(u_bound - l_bound >= n);

        if (n == 1) return { randomInt(l_bound, u_bound - 1) };

        std::vector<IntType> nums(u_bound - l_bound);
        std::iota(nums.begin(), nums.end(), l_bound);  // [l_bound, u_bound)

        for (size_t i = 0; i < n; i++)
        {
            size_t idx = randomInt(i, nums.size() - 1);
            std::swap(nums[idx], nums[i]);
        }
        nums.resize(n);

        return nums;
    }

    std::vector<double> randomSimplexPoint(size_t dim)
    {
        assert(dim > 0);

        std::uniform_real_distribution dist{ 0.0, 1.0 };

        std::vector<double> point;
        point.reserve(dim);

        double sum = 0.0;
        for (size_t i = 0; i < dim; i++)
        {
            point.push_back(-std::log(dist(prng)));
            sum += point.back();
        }
        for (size_t i = 0; i < dim; i++)
        {
            point[i] /= sum;
        }

        return point;
    }

    size_t sampleCdf(const std::vector<double>& cdf)
    {
        assert(!cdf.empty());

        auto it = std::lower_bound(cdf.begin(), cdf.end(), randomReal(0.0, cdf.back())); // cdf.back() in case not exactly 1.0

        return size_t(it - cdf.begin());
    }

}

#endif // !GA_RANDOM_HPP