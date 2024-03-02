/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_RNG_HPP
#define GA_UTILITY_RNG_HPP

#include "small_vector.hpp"
#include "type_traits.hpp"
#include "bit.hpp"
#include "rcu.hpp"
#include "indestructible.hpp"
#include "utility.hpp"
#include <algorithm>
#include <functional>
#include <ranges>
#include <array>
#include <vector>
#include <span>
#include <unordered_set>
#include <bit>
#include <random>
#include <limits>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <cstdint>
#include <cstddef>
#include <concepts>


#ifndef GAPP_SEED
#   define GAPP_SEED 0x3da99432ab975d26LL
#endif


/** Contains the PRNG classes and functions used for generating random numbers. */
namespace gapp::rng
{
    /**
    * Splitmix64 pseudo-random number generator based on
    * https://prng.di.unimi.it/splitmix64.c. This generator
    * is only used for seeding the Xoroshiro generators.
    */
    class Splitmix64
    {
    public:
        using result_type = std::uint64_t;  /**< The generator generates 64 bit integers. */
        using state_type  = std::uint64_t;  /**< The generator has a 64 bit state. */

        /**
        * Create a new Splitmix64 generator initialized from a 64 bit seed value.
        * 
        * @note
        *  The global prng instance should be used instead of creating
        *  new instances of this generator.
        * 
        * @param seed The seed used to initialize the state of the generator.
        */
        explicit constexpr Splitmix64(state_type seed) noexcept :
            state_(seed)
        {}

        /** @returns The next number of the sequence. */
        constexpr result_type operator()() noexcept
        {
            result_type z = (state_ += 0x9e3779b97f4a7c15);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;

            return z ^ (z >> 31);
        }

        /** Set a new seed for the generator. */
        constexpr void seed(state_type seed) noexcept { state_ = seed; }

        /** @returns The smallest possible value that can be generated. */
        static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }

        /** @returns The largest possible value that can be generated. */
        static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

        /** Compare the internal state of 2 generators. @returns true if they are the same. */
        friend constexpr bool operator==(const Splitmix64&, const Splitmix64&) = default;

    private:
        state_type state_;
    };


    /**
    * Xoroshiro128+ pseudo-random number generator based on
    * https://prng.di.unimi.it/xoroshiro128plus.c
    * 
    * @see
    *   David Blackman and Sebastiano Vigna. "Scrambled linear pseudorandom number generators."
    *   ACM Transactions on Mathematical Software 47, no. 4 (2021): 1-32.
    */
    class Xoroshiro128p
    {
    public:
        using result_type = std::uint64_t;           /**< The generator generates 64 bit integers. */
        using state_type  = std::array<uint64_t, 2>; /**< The generator has 128 bits of state. */

        /**
        * Create a Xoroshiro128+ generator initialized from a 64 bit seed
        * value.
        *
        * @note
        *  The global prng instance should be used instead of creating
        *  new instances of this generator.
        *
        * @param seed The seed used to initialize the state of the generator.
        */
        explicit constexpr Xoroshiro128p(std::uint64_t seed) noexcept :
            state_(seed_sequence(seed))
        {}

        /** @returns The next number of the sequence. */
        constexpr result_type operator()() noexcept
        {
            const auto result = state_[0] + state_[1];
            const auto xstate = state_[0] ^ state_[1];

            state_[0] = std::rotl(state_[0], 24) ^ xstate ^ (xstate << 16);
            state_[1] = std::rotl(xstate, 37);

            return result;
        }

        /** Advance the state of the generator by 2^96 steps. */
        constexpr Xoroshiro128p& jump() noexcept
        {
            state_type new_state{ 0, 0 };

            for (std::uint64_t JUMP : { 0xd2a98b26625eee7bULL, 0xdddf9b1090aa7ac1ULL })
            {
                for (std::size_t n = 0; n < 64; n++)
                {
                    if (detail::is_nth_bit_set(JUMP, n))
                    {
                        new_state[0] ^= state_[0];
                        new_state[1] ^= state_[1];
                    }
                    std::invoke(*this);
                }
            }
            state_ = new_state;

            return *this;
        }

        /** Set a new seed for the generator. */
        constexpr void seed(std::uint64_t seed) noexcept { state_ = seed_sequence(seed); }

        /** @returns The smallest possible value that can be generated. */
        static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }

        /** @returns The largest possible value that can be generated. */
        static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

        /** Compare the internal state of 2 generators. @returns True if they are the same. */
        friend constexpr bool operator==(const Xoroshiro128p&, const Xoroshiro128p&) = default;

    private:
        static constexpr state_type seed_sequence(std::uint64_t seed) noexcept
        {
            Splitmix64 seed_seq_gen(seed);
            return { seed_seq_gen(), seed_seq_gen() };
        }

        alignas(128) state_type state_;
    };


    /**
     * The pseudo-random number generator class used in the library.
     * This class is a simple wrapper around the Xoroshiro128p generator
     * to make it thread-safe.
     * 
     * @note
     *  The global prng instance should be used instead of creating
     *  new instances of this generator.
     */
    class ConcurrentXoroshiro128p
    {
    public:
        using result_type = Xoroshiro128p::result_type;
        using state_type  = Xoroshiro128p::state_type;

        /** @return The next number of the sequence. Thread-safe. */
        result_type operator()() const noexcept
        {
            //std::shared_lock _{ generator_.instance };
            return std::invoke(*generator_.instance);
        }

        /** Set a new seed for the generator. Thread-safe. */
        static void seed(std::uint64_t seed)
        {
            std::scoped_lock _{ tls_generators_->lock };
            global_generator_.seed(seed);
            for (Generator* generator : tls_generators_->list)
            {
                *generator = global_generator_.jump();
            }
        }

        /** @returns The smallest possible value that can be generated. */
        static constexpr result_type min() noexcept { return Xoroshiro128p::min(); }

        /** @returns The largest possible value that can be generated. */
        static constexpr result_type max() noexcept { return Xoroshiro128p::max(); }

    private:
        using Generator = detail::rcu_obj<Xoroshiro128p>;

        struct RegisteredGenerator
        {
            RegisteredGenerator() noexcept // NOLINT(*exception-escape)
            {
                std::scoped_lock _{ tls_generators_->lock };
                instance.get() = global_generator_.jump();
                tls_generators_->list.push_back(std::addressof(instance));
            }

            ~RegisteredGenerator() noexcept
            {
                std::scoped_lock _{ tls_generators_->lock };
                std::erase(tls_generators_->list, std::addressof(instance));
            }

            Generator instance{ 0 };
        };

        struct GeneratorList
        {
            detail::spinlock lock;
            std::vector<Generator*> list;
        };

        GAPP_API inline static Xoroshiro128p global_generator_{ GAPP_SEED };
        GAPP_API inline static detail::Indestructible<GeneratorList> tls_generators_;
        alignas(128) inline static thread_local RegisteredGenerator generator_;
    };


    /** The global pseudo-random number generator instance used in the algorithms. */
    inline constexpr ConcurrentXoroshiro128p prng;


    /** Generate a random boolean value from a uniform distribution. */
    inline bool randomBool() noexcept;

    /** Generate a random integer from a uniform distribution on the closed interval [lbound, ubound]. */
    template<std::integral IntType = int>
    IntType randomInt(IntType lbound, IntType ubound);

    /** Generate a random floating-point value from a uniform distribution on the half-open interval [0.0, 1.0). */
    template<std::floating_point RealType = double>
    RealType randomReal();

    /** Generate a random floating-point value of from a uniform distribution on the half-open interval [lbound, ubound). */
    template<std::floating_point RealType = double>
    RealType randomReal(RealType lbound, RealType ubound);

    /** Generate a random floating-point value from a normal distribution with the specified mean and std deviation. */
    template<std::floating_point RealType = double>
    RealType randomNormal(RealType mean = 0.0, RealType std_dev = 1.0);

    /** Generate a random integer value from an approximate binomial distribution with the parameters n and p. */
    template<std::integral IntType = int>
    IntType randomBinomial(IntType n, double p);


    /** Generate a random index for a range. */
    template<std::ranges::random_access_range R>
    detail::size_type<R> randomIndex(const R& range);

    /** Pick a random element from a range. */
    template<std::forward_iterator Iter>
    Iter randomElement(Iter first, Iter last);


    /** Generate @p count unique integers from the half-open range [@p lbound, @p ubound). */
    template<std::integral IntType>
    small_vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count);


    /** Select an index based on a discrete CDF. */
    template<std::floating_point T>
    size_t sampleCdf(std::span<const T> cdf);

} // namespace gapp::rng


/* IMPLEMENTATION */

#include <utility>
#include <numeric>
#include <iterator>
#include <cmath>
#include <climits>

namespace gapp::rng
{
    bool randomBool() noexcept
    {
        constinit thread_local uint64_t bit_pool = 1;

        if (bit_pool == detail::lsb_mask<uint64_t>)
        {
            bit_pool = prng() | detail::msb_mask<uint64_t>;
        }

        const bool bit = detail::last_bit(bit_pool);
        bit_pool >>= 1;

        return bit;
    }

    template<std::integral IntType>
    IntType randomInt(IntType lbound, IntType ubound)
    {
        GAPP_ASSERT(lbound <= ubound);

        // std::uniform_int_distribution doesnt support char and other 8 bit types
        using NonCharInt = std::conditional_t<std::is_signed_v<IntType>, int64_t, uint64_t>;

        std::uniform_int_distribution<NonCharInt> dist{ lbound, ubound };

        return dist(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal()
    {
        return (rng::prng() >> detail::exponent_bits<RealType>) * (std::numeric_limits<RealType>::epsilon() / 2.0);
    }

    template<std::floating_point RealType>
    RealType randomReal(RealType lbound, RealType ubound)
    {
        GAPP_ASSERT(lbound <= ubound);

        return lbound + (ubound - lbound) * rng::randomReal();
    }

    template<std::floating_point RealType>
    RealType randomNormal(RealType mean, RealType std_dev)
    {
        GAPP_ASSERT(std_dev >= 0.0);

        // keep the distribution for the state
        thread_local std::normal_distribution<RealType> dist;

        return std_dev * dist(rng::prng) + mean;
    }

    template<std::integral IntType>
    IntType randomBinomialApproxNormal(IntType n, double p)
    {
        GAPP_ASSERT(n >= 0);
        GAPP_ASSERT(0.0 <= p && p <= 1.0);

        const double mean = n * p;
        const double std_dev = std::sqrt(mean * (1.0 - p));

        const double accept_min = -0.5;
        const double accept_max = n + 0.5;

        double rand = rng::randomNormal(mean, std_dev);
        while (!(accept_min < rand && rand < accept_max))
        {
            rand = rng::randomNormal(mean, std_dev);
        }

        return static_cast<IntType>(std::round(rand));
    }

    template<std::integral IntType>
    IntType randomBinomialApproxPoisson(IntType n, double p)
    {
        GAPP_ASSERT(n >= 0);
        GAPP_ASSERT(0.0 <= p && p <= 1.0);
        GAPP_ASSERT(n * p < 12.0); // required to avoid data race in std::poisson_distribution with libstdc++
        
        if (p == 0.0) return 0;
        if (p == 1.0) return n;

        // keep the distribution to avoid expensive initialization
        thread_local std::poisson_distribution<IntType> dist;

        if (dist.mean() != n * p)
        {
            dist = std::poisson_distribution<IntType>{ n * p };
        }

        IntType rand = dist(rng::prng);
        while (rand > n) { rand = dist(rng::prng); }
        return rand;
    }

    template<std::integral IntType>
    IntType randomBinomialExact(IntType n, double p)
    {
        GAPP_ASSERT(n >= 0);
        GAPP_ASSERT(0.0 <= p && p <= 1.0);

        IntType res = 0;
        while (n--) { res += (randomReal() <= p); }
        return res;
    }

    template<std::integral IntType>
    IntType randomBinomial(IntType n, double p)
    {
        GAPP_ASSERT(n >= 0);
        GAPP_ASSERT(0.0 <= p && p <= 1.0);

        const double mean = p * n;

        const bool use_poisson = (n > 50 && mean < 10.0) || ((5 < n && n <= 50) && p < 0.2);
        const bool use_normal = !use_poisson && (n > 5);

        if (use_poisson) return rng::randomBinomialApproxPoisson(n, p);
        if (use_normal)  return rng::randomBinomialApproxNormal(n, p);
        return rng::randomBinomialExact(n, p);
    }

    template<std::ranges::random_access_range R>
    detail::size_type<R> randomIndex(const R& range)
    {
        GAPP_ASSERT(!range.empty());

        return std::uniform_int_distribution<detail::size_type<R>>{ 0, range.size() - 1 }(rng::prng);
    }

    template<std::forward_iterator Iter>
    Iter randomElement(Iter first, Iter last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);

        const auto max_offset = std::distance(first, last) - 1;
        const auto offset = rng::randomInt(0_pd, max_offset);

        return std::next(first, offset);
    }

    template<std::integral IntType>
    GAPP_NOINLINE small_vector<IntType> sampleUniqueSet(IntType lbound, IntType ubound, size_t count)
    {
        std::unordered_set<IntType> selected(count);
        small_vector<IntType> numbers(count);

        IntType limit = ubound - detail::promoted_t<IntType>(count);

        for (auto number = numbers.begin(); limit < ubound; ++limit, ++number)
        {
            const auto n = rng::randomInt(lbound, limit);
            *number = selected.contains(n) ? limit : n;
            selected.insert(n);
        }

        return numbers;
    }

    template<std::integral IntType>
    small_vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count)
    {
        const size_t range_len = detail::range_length(lbound, ubound);

        GAPP_ASSERT(ubound >= lbound);
        GAPP_ASSERT(range_len >= count);

        const bool select_many = (count > 0.6 * range_len);
        const bool huge_range  = (range_len >= 65536);

        if (huge_range) [[unlikely]] return rng::sampleUniqueSet(lbound, ubound, count);

        small_vector<IntType> numbers(count);

        thread_local std::vector<bool> is_selected;
        is_selected.resize(range_len);
        std::fill(is_selected.begin(), is_selected.end(), select_many);

        if (!select_many)
        {
            IntType limit = ubound - detail::promoted_t<IntType>(count);

            for (auto number = numbers.begin(); limit < ubound; ++limit, ++number)
            {
                const auto n = rng::randomInt(lbound, limit);
                *number = is_selected[n - lbound] ? limit : n;
                is_selected[*number - lbound] = true;
            }
        }
        else
        {
            IntType rcount = detail::promoted_t<IntType>(range_len - count);

            for (IntType limit = ubound - rcount; limit < ubound; ++limit)
            {
                const auto n = rng::randomInt(lbound, limit);
                (is_selected[n - lbound] ? is_selected[n - lbound] : is_selected[limit - lbound]) = false;
            }

            auto number = numbers.begin();
            for (size_t n = 0; n < is_selected.size(); n++)
            {
                if (is_selected[n]) { *number++ = lbound + n; }
            }
        }

        return numbers;
    }

    template<std::floating_point T>
    size_t sampleCdf(std::span<const T> cdf)
    {
        GAPP_ASSERT(!cdf.empty());
        GAPP_ASSERT(0.0 <= cdf.front());

        const auto limitval = rng::randomReal<T>(0.0, cdf.back()); // use cdf.back() in case it's not exactly 1.0
        const auto selected = std::lower_bound(cdf.begin(), cdf.end(), limitval);

        return static_cast<size_t>(selected - cdf.begin());
    }

} // namespace gapp::rng

#endif // !GA_UTILITY_RNG_HPP