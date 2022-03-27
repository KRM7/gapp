/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_RANDOM_HPP
#define GA_RANDOM_HPP

#include "utils.hpp"
#include "concepts.hpp"

#include <vector>
#include <random>
#include <cstdint>
#include <cstddef>
#include <concepts>

/** Contains the PRNG classes and functions for generating random numbers. */
namespace genetic_algorithm::rng
{
    /**
    * Splitmix64 PRNG adapted from https://prng.di.unimi.it/splitmix64.c \n
    * Only used for seeding the other PRNGs.
    */
    class Splitmix64
    {
    public:
        using result_type = uint_fast64_t;
        using state_type = uint_fast64_t;

        explicit Splitmix64(state_type seed);

        result_type operator()() noexcept;

    private:
        state_type state_;
    };

    /** xoroshiro128+ PRNG adapted from https://prng.di.unimi.it/xoroshiro128plus.c */
    class Xoroshiro128p
    {
    public:
        using result_type = uint_fast64_t;
        using state_type = uint_fast64_t;

        explicit Xoroshiro128p(uint_fast64_t seed);

        result_type operator()() noexcept;

        static constexpr result_type min() noexcept;
        static constexpr result_type max() noexcept;

    private:
        state_type state_[2];

        static state_type rotl64(state_type value, unsigned shift) noexcept;
    };

    /** Thread-safe seed generator for seeding PRNGs created on different threads. */
    class SeedGenerator
    {
    public:
        /** Generate a new seed that can be used to initialize a PRNG. */
        Splitmix64::result_type operator()();
    private:
        Splitmix64 gen_ = Splitmix64{ GA_SEED() };
    };

    /** Global seed generator used in the genetic algorithms to seed PRNGs. */
    inline SeedGenerator seed_gen{};

    /** The PRNG type used in the genetic algorithms. */
    using PRNG = Xoroshiro128p;

    /** PRNG instance(s) used in the genetic algorithms. */
    thread_local inline PRNG prng{ seed_gen() };

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

    /** Generates a random integer of type IntType from a uniform distribution on the closed interval [l_bound, u_bound]. */
    template<std::integral IntType = int>
    IntType randomInt(IntType l_bound, IntType u_bound);

    /** Generates a random idx for a container. */
    template<detail::IndexableContainer T>
    size_t randomIdx(const T& cont);

    /** Pick a random element from a container. */
    template<detail::Container T>
    auto randomElement(const T& cont) -> typename T::value_type;

    /** Generates a random boolean value from a uniform distribution. */
    inline bool randomBool() noexcept;

    /** Generates n unique integers from the range [0, u_bound). */
    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType u_bound, size_t n);

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
    inline Splitmix64::Splitmix64(state_type seed)
        : state_(seed)
    {
    }

    inline Splitmix64::result_type Splitmix64::operator()() noexcept
    {
        state_ += 0x9e3779b97f4a7c15;
        result_type z = state_;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;

        return z ^ (z >> 31);
    }


    inline Xoroshiro128p::Xoroshiro128p(uint_fast64_t seed)
    {
        Splitmix64 seed_seq_gen(seed);
        state_[0] = seed_seq_gen();
        state_[1] = seed_seq_gen();
    }

    inline Xoroshiro128p::result_type Xoroshiro128p::operator()() noexcept
    {
        state_type s0 = state_[0];
        state_type s1 = state_[1];
        result_type result = s0 + s1;

        s1 ^= s0;
        state_[0] = rotl64(s0, 24) ^ s1 ^ (s1 << 16);
        state_[1] = rotl64(s1, 37);

        return result;
    }

    inline constexpr Xoroshiro128p::result_type Xoroshiro128p::min() noexcept
    {
        return std::numeric_limits<result_type>::lowest();
    }

    inline constexpr Xoroshiro128p::result_type Xoroshiro128p::max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

    inline Xoroshiro128p::state_type Xoroshiro128p::rotl64(state_type value, unsigned shift) noexcept
    {
        return (value << shift) | (value >> (64U - shift));
    }


    inline Splitmix64::result_type SeedGenerator::operator()()
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lock(m);

        return gen_();
    }


    template<std::floating_point RealType>
    RealType randomReal()
    {
        std::uniform_real_distribution<RealType> distribution{ 0.0, 1.0 };

        return distribution(prng);
    }

    template<std::floating_point RealType>
    RealType randomReal(RealType l_bound, RealType u_bound)
    {
        assert(l_bound <= u_bound);

        std::uniform_real_distribution<RealType> distribution{ l_bound, u_bound };

        return distribution(prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal()
    {
        std::normal_distribution<RealType> distribution{ 0.0, 1.0 };

        return distribution(prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal(RealType mean, RealType SD)
    {
        assert(SD > 0.0);

        std::normal_distribution<RealType> distribution{ mean, SD };

        return distribution(prng);
    }

    template<std::integral IntType>
    IntType randomInt(IntType l_bound, IntType u_bound)
    {
        assert(l_bound <= u_bound);

        std::uniform_int_distribution<IntType> distribution{ l_bound, u_bound };

        return distribution(prng);
    }

    template<detail::IndexableContainer T>
    size_t randomIdx(const T& cont)
    {
        assert(!cont.empty()); /* There are no valid indices otherwise. */

        std::uniform_int_distribution<size_t> distribution{ 0, cont.size() - 1 };

        return distribution(prng);
    }

    template<detail::Container T>
    auto randomElement(const T& cont) -> typename T::value_type
    {
        size_t i = randomInt( size_t{0}, cont.size() - 1 );

        return *std::next(cont.begin(), i);
    }

    bool randomBool() noexcept
    {
        std::uniform_int_distribution distribution{ 0, 1 };

        return distribution(prng);
    }

    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType u_bound, size_t n)
    {
        assert(u_bound >= n);

        std::vector<IntType> nums(u_bound);
        std::iota(nums.begin(), nums.end(), IntType{ 0 });  // [0, u_bound)

        for (size_t i = 0; i < n; i++)
        {
            size_t idx = randomInt(i, nums.size() - 1);
            std::swap(nums[idx], nums[i]);
        }
        nums.resize(n);

        return nums;
    }

    template<std::integral IntType>
    std::vector<IntType> sampleUnique(IntType l_bound, IntType u_bound, size_t n)
    {
        assert(l_bound <= u_bound);
        assert(u_bound - l_bound >= n);

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

        static thread_local std::minstd_rand0 rng{ GA_SEED() };
        std::uniform_real_distribution<double> dist{ 0.0, 1.0 };

        std::vector<double> point;
        point.reserve(dim);

        double sum = 0.0;
        for (size_t i = 0; i < dim; i++)
        {
            point.push_back(-std::log(dist(rng)));
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

        auto idx = std::lower_bound(cdf.begin(), cdf.end(), rng::randomReal(0.0, cdf.back())) - cdf.begin();

        return size_t(idx);
    }

}

#endif // !GA_RANDOM_HPP