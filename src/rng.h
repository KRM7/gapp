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
* in the genetic algorithms.
*
* @file rng.h
*/

#ifndef GA_RANDOM_H
#define GA_RANDOM_H

#include <random>
#include <limits>
#include <vector>
#include <cmath>
#include <cstdint>
#include <cassert>


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

        splitmix64(state_type seed) : state(seed) {}

        result_type operator()() noexcept
        {
            state += 0x9e3779b97f4a7c15;
            result_type z = state;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            return z ^ (z >> 31);
        }

    private:
        state_type state;
    };

    /** xoroshiro128+ PRNG adapted from https://prng.di.unimi.it/xoroshiro128plus.c */
    class xoroshiro128p
    {
    public:
        using result_type = uint_fast64_t;
        using state_type = uint_fast64_t;

        xoroshiro128p(uint_fast64_t seed)
        {
            splitmix64 seed_seq_gen(seed);
            state[0] = seed_seq_gen();
            state[1] = seed_seq_gen();
        }

        result_type operator()() noexcept
        {
            state_type s0 = state[0];
            state_type s1 = state[1];
            result_type result = s0 + s1;

            s1 ^= s0;
            state[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16);
            state[1] = rotl(s1, 37);

            return result;
        }

        static constexpr result_type min() noexcept
        {
            return std::numeric_limits<result_type>::lowest();
        }
        static constexpr result_type max() noexcept
        {
            return std::numeric_limits<result_type>::max();
        }

    private:
        state_type state[2];

        static state_type rotl(state_type x, int k) noexcept
        {
            return (x << k) | (x >> (64 - k));
        }
    };

    /** The PRNG used in the genetic algorithm. */
    using PRNG = xoroshiro128p;

    /** Generates a random double on the interval [l_bound, u_bound). */
    inline double generateRandomDouble(double l_bound = 0.0, double u_bound = 1.0)
    {
        assert(l_bound <= u_bound);

        static thread_local PRNG engine{ std::random_device{}() };
        std::uniform_real_distribution<double> distribution{ l_bound, u_bound };

        return distribution(engine);
    }

    /** Generates a random integer of type T on the closed interval [l_bound, u_bound]. */
    template<typename T>
    inline T generateRandomInt(T l_bound, T u_bound)
    {
        assert(l_bound <= u_bound);

        static thread_local PRNG engine{ std::random_device{}() };
        std::uniform_int_distribution<T> distribution{ l_bound, u_bound };

        return distribution(engine);
    }

    /**
    * Generates a random unsigned integer on the closed interval [0, c_size-1]. \n
    * Used to generate a random index for containers, with c_size being the size of the container.
    */
    inline size_t generateRandomIdx(size_t c_size)
    {
        return generateRandomInt(size_t{ 0 }, c_size - 1);
    }

    /** Generates a random boolean value. */
    inline bool generateRandomBool()
    {
        return bool(generateRandomInt(size_t{ 0 }, size_t{ 1 }));
    }

    /** Generates a random double from a normal distribution with the parameters mean and SD. */
    inline double generateRandomNorm(double mean = 0.0, double SD = 1.0)
    {
        assert(SD > 0.0);

        static thread_local PRNG engine{ std::random_device{}() };
        std::normal_distribution<double> distribution{ mean, SD };

        return distribution(engine);
    }

    /**
    * Sample a point from a uniform distribution on a unit simplex in @p dim dimensions. \n
    *
    * @param dim The number of elements in the generated vector.
    * @returns The generated random point/vector.
    */
    inline std::vector<double> generateRandomSimplexPoint(size_t dim)
    {
        assert(dim > 0);

        static thread_local std::minstd_rand0 engine{ std::random_device{}() };
        std::uniform_real_distribution<double> distribution{ 0.0, 1.0 };

        std::vector<double> vec;
        vec.reserve(dim);

        double sum = 0.0;
        for (size_t i = 0; i < dim; i++)
        {
            vec.push_back(-std::log(distribution(engine)));
            sum += vec.back();
        }
        for (size_t i = 0; i < dim; i++)
        {
            vec[i] /= sum;
        }

        return vec;
    }

} // namespace genetic_algorithm::rng

#endif // !GA_RANDOM_H