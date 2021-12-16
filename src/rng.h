/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

/**
* This file contains the PRNG classes and functions used for generating random numbers
* in the genetic algorithms. The rng functions are thread safe.
*
* @file rng.h
*/

#ifndef GA_RANDOM_H
#define GA_RANDOM_H

#include <vector>
#include <random>
#include <cstdint>
#include <cstddef>

/** Contains the PRNG classes and functions for generating random numbers. */
namespace genetic_algorithm::rng
{
    /**
    * Splitmix64 PRNG adapted from https://prng.di.unimi.it/splitmix64.c \n
    * Only used for seeding the other PRNGs.
    */
    class splitmix64
    {
    public:
        using result_type = uint_fast64_t;
        using state_type = uint_fast64_t;

        explicit splitmix64(state_type seed);

        result_type operator()() noexcept;

    private:
        state_type state;
    };

    /** xoroshiro128+ PRNG adapted from https://prng.di.unimi.it/xoroshiro128plus.c */
    class xoroshiro128p
    {
    public:
        using result_type = uint_fast64_t;
        using state_type = uint_fast64_t;

        explicit xoroshiro128p(uint_fast64_t seed);

        result_type operator()() noexcept;

        static constexpr result_type min() noexcept;
        static constexpr result_type max() noexcept;

    private:
        state_type state[2];

        static state_type rotl(state_type x, int k) noexcept;
    };

    /** Thread-safe seed generator for seeding PRNGs created on different threads. */
    class SeedGenerator
    {
    public:
        /** Generate a new seed that can be used to initialize a PRNG. */
        splitmix64::result_type operator()() noexcept;
    private:
        splitmix64 gen_ = splitmix64{ ::std::random_device{}() };
    };

    /** Global seed generator used in the genetic algorithms to seed PRNGs. */
    inline SeedGenerator seed_gen{};

    /** The PRNG type used in the genetic algorithms. */
    using PRNG = xoroshiro128p;

    /** PRNG instance(s) used in the genetic algorithm. */
    thread_local inline PRNG prng{ seed_gen() };

    /** Generates a random floating-point value of type RealType from a uniform distribution on the interval [0.0, 1.0). */
    template<typename RealType = double>
    inline RealType randomReal();

    /** Generates a random floating-point value of type RealType from a uniform distribution on the interval [l_bound, u_bound). */
    template<typename RealType = double>
    inline RealType randomReal(RealType l_bound, RealType u_bound);

    /** Generates a random floating-point value of type RealType from a standard normal distribution. */
    template<typename RealType = double>
    inline RealType randomNormal();

    /** Generates a random floating-point value of type RealType from a normal distribution with the parameters mean and SD. */
    template<typename RealType = double>
    inline RealType randomNormal(RealType mean, RealType SD);

    /** Generates a random integer of type IntType from a uniform distribution on the closed interval [l_bound, u_bound]. */
    template<typename IntType = int>
    inline IntType randomInt(IntType l_bound, IntType u_bound);

    /**
    * Generates a random unsigned integer from a uniform distribution on the closed interval [0, c_size-1]. \n
    * Used to generate a random index for containers, with c_size being the size of the container.
    */
    inline size_t randomIdx(size_t c_size);

    /** Generates a random boolean value from a uniform distribution. */
    inline bool randomBool();

    /** Generates n unique integers from the range [0, u_bound). */
    template<typename IntType>
    inline std::vector<IntType> sampleUnique(IntType u_bound, size_t n);

} // namespace genetic_algorithm::rng


/* IMPLEMENTATION */

#include <limits>
#include <utility>
#include <numeric>
#include <mutex>
#include <cassert>

namespace genetic_algorithm::rng
{
    inline splitmix64::splitmix64(state_type seed)
        : state(seed)
    {
    }

    inline splitmix64::result_type splitmix64::operator()() noexcept
    {
        state += 0x9e3779b97f4a7c15;
        result_type z = state;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;

        return z ^ (z >> 31);
    }


    inline xoroshiro128p::xoroshiro128p(uint_fast64_t seed)
    {
        splitmix64 seed_seq_gen(seed);
        state[0] = seed_seq_gen();
        state[1] = seed_seq_gen();
    }

    inline xoroshiro128p::result_type xoroshiro128p::operator()() noexcept
    {
        state_type s0 = state[0];
        state_type s1 = state[1];
        result_type result = s0 + s1;

        s1 ^= s0;
        state[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16);
        state[1] = rotl(s1, 37);

        return result;
    }

    inline constexpr xoroshiro128p::result_type xoroshiro128p::min() noexcept
    {
        return std::numeric_limits<result_type>::lowest();
    }

    inline constexpr xoroshiro128p::result_type xoroshiro128p::max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }

    inline xoroshiro128p::state_type xoroshiro128p::rotl(state_type x, int k) noexcept
    {
        return (x << k) | (x >> (64 - k));
    }

    inline splitmix64::result_type SeedGenerator::operator()()
    inline splitmix64::result_type SeedGenerator::operator()() noexcept
    {
        static std::mutex m;
        std::lock_guard<std::mutex> lock(m);

        return gen_();
    }


    template<typename RealType>
    inline RealType randomReal()
    {
        std::uniform_real_distribution<RealType> distribution{ 0.0, 1.0 };

        return distribution(prng);
    }

    template<typename RealType>
    RealType randomReal(RealType l_bound, RealType u_bound)
    {
        assert(l_bound <= u_bound);

        std::uniform_real_distribution<RealType> distribution{ l_bound, u_bound };

        return distribution(prng);
    }

    template<typename RealType>
    RealType randomNormal()
    {
        std::normal_distribution<RealType> distribution{ 0.0, 1.0 };

        return distribution(prng);
    }

    template<typename RealType>
    RealType randomNormal(RealType mean, RealType SD)
    {
        assert(SD > 0.0);

        std::normal_distribution<RealType> distribution{ mean, SD };

        return distribution(prng);
    }

    template<typename IntType>
    IntType randomInt(IntType l_bound, IntType u_bound)
    {
        assert(l_bound <= u_bound);

        std::uniform_int_distribution<IntType> distribution{ l_bound, u_bound };

        return distribution(prng);
    }

    size_t randomIdx(size_t c_size)
    {
        assert(c_size > 0); /* There are no valid indices otherwise. */

        std::uniform_int_distribution<size_t> distribution{ 0, c_size - 1 };

        return distribution(prng);
    }

    bool randomBool()
    {
        std::uniform_int_distribution<int> distribution{ 0, 1 };

        return distribution(prng);
    }

    template<typename IntType>
    std::vector<IntType> sampleUnique(IntType u_bound, size_t n)
    {
        std::vector<IntType> nums(u_bound);
        std::iota(nums.begin(), nums.end(), IntType{ 0 });

        for (size_t i = 0; i < n; i++)
        {
            size_t idx = randomIdx(nums.size() - i);
            std::swap(nums[idx], *(nums.end() - 1 - i));
        }

        return std::vector(nums.end() - n, nums.end());
    }

}

#endif // !GA_RANDOM_H