/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_RNG_HPP
#define GAPP_UTILITY_RNG_HPP

#include "algorithm.hpp"
#include "functional.hpp"
#include "thread_pool.hpp"
#include "distribution.hpp"
#include "small_vector.hpp"
#include "dynamic_bitset.hpp"
#include "type_traits.hpp"
#include "bit.hpp"
#include "spinlock.hpp"
#include "utility.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <ranges>
#include <array>
#include <vector>
#include <span>
#include <unordered_set>
#include <concepts>
#include <bit>
#include <random>
#include <limits>
#include <memory>
#include <mutex>
#include <utility>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <climits>


#ifndef GAPP_SEED
#   define GAPP_SEED 0x3da99432ab975d26LL
#endif


/** Contains the PRNG classes and functions used for generating random numbers. */
namespace gapp::rng
{
    /** Generate a random boolean value from a uniform distribution. Thread-safe. */
    inline bool randomBool() noexcept;

    /** Generate a random integer from a uniform distribution on the closed interval [lbound, ubound]. Thread-safe. */
    template<std::integral IntType = int>
    IntType randomInt(IntType lbound, IntType ubound) noexcept;

    /** Generate a random floating-point value from a uniform distribution on the half-open interval [0.0, 1.0). Thread-safe. */
    template<std::floating_point RealType = double>
    RealType randomReal() noexcept;

    /** Generate a random floating-point value of from a uniform distribution on the half-open interval [lbound, ubound). Thread-safe. */
    template<std::floating_point RealType = double>
    RealType randomReal(RealType lbound, RealType ubound) noexcept;

    /** Generate a random floating-point value from a normal distribution with the specified mean and std deviation. Thread-safe. */
    template<std::floating_point RealType = double>
    RealType randomNormal(RealType mean = 0.0, RealType std_dev = 1.0) noexcept;

    /** Generate a random integer value from a binomial distribution with the parameters n and p. Thread-safe. */
    template<std::integral IntType = int>
    IntType randomBinomial(IntType n, double p) noexcept;


    /** Generate random integer values from a binomial distribution. operator() is thread-safe. */
    template<std::integral IntType = int>
    class CachedRandomBinomial
    {
    public:
        void init(IntType n, double p) noexcept;
        IntType operator()(IntType n, double p) const noexcept;
    private:
        mutable rng::binomial_distribution<IntType> dist_ = { 0, 0.0 };
    };


    /** Generate a random index for a range (from a uniform distribution over the valid indices). Thread-safe. */
    template<std::ranges::random_access_range R>
    detail::size_type<R> randomIndex(const R& range) noexcept;

    /** Pick a random element from a range (from a uniform distribution over the elements). Thread-safe. */
    template<std::forward_iterator Iter>
    Iter randomElement(Iter first, Iter last) noexcept;

    /** Pick a random element from a range using the specified cumulative distribution function. Thread-safe. */
    template<std::ranges::random_access_range Range, std::ranges::random_access_range Dist>
    std::ranges::range_reference_t<Range> randomElement(Range&& range, Dist&& cdf);


    /** Generate @p count unique integers from the half-open range [@p lbound, @p ubound). Thread-safe. */
    template<std::integral IntType>
    small_vector<IntType> sampleUnique(IntType lbound, IntType ubound, size_t count);

    /** Select an index based on a discrete CDF. */
    template<std::ranges::random_access_range Range>
    size_t sampleCdf(const Range& cdf);


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

        state_type state_;
    };


    /**
     * The pseudo-random number generator class used in the library.
     * This class is a simple wrapper around the Xoroshiro128p generator
     * to make operator() thread-safe.
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
        GAPP_API result_type operator()() const noexcept
        {
            return std::invoke(generator_.instance);
        }

        /** 
         * Set a new seed for the generator. This function is not thread-safe and shouldn't
         * be called concurrently with the random number generation functions (e.g. while a GA
         * is running).
         */
        GAPP_API static void seed(std::uint64_t seed)
        {
            detail::parallel_for(detail::iota_iterator(0_sz), detail::iota_iterator(execution_threads()), [](size_t)
            {
                std::ignore = generator_;
            });

            std::scoped_lock _{ tls_generators()->lock};
            global_generator().seed(seed);

            std::ranges::sort(tls_generators()->list, std::greater{}, [](RegisteredGenerator* generator)
            {
                return generator->thread_id->load(std::memory_order_acquire);
            });

            for (RegisteredGenerator* generator : tls_generators()->list)
            {
                generator->instance = global_generator().jump();
            }
        }

        /** @returns The smallest possible value that can be generated. */
        static constexpr result_type min() noexcept { return Xoroshiro128p::min(); }

        /** @returns The largest possible value that can be generated. */
        static constexpr result_type max() noexcept { return Xoroshiro128p::max(); }

        GAPP_API static auto& generator() noexcept { return generator_; }

    private:
        struct RegisteredGenerator
        {
            RegisteredGenerator() noexcept // NOLINT(*exception-escape)
            {
                std::scoped_lock _{ tls_generators()->lock };
                instance = global_generator().jump();
                thread_id = &detail::thread_pool::this_thread_id();
                tls_generators()->list.push_back(this);
            }

            ~RegisteredGenerator() noexcept
            {
                std::scoped_lock _{ tls_generators()->lock };
                std::erase(tls_generators()->list, this);
            }

            Xoroshiro128p instance{ 0 };
            std::atomic<std::uint64_t>* thread_id;
        };

        struct GeneratorList
        {
            detail::spinlock lock;
            std::vector<RegisteredGenerator*> list;
        };

        GAPP_API static Xoroshiro128p& global_generator() noexcept
        {
            static Xoroshiro128p global_generator_{ GAPP_SEED };
            return global_generator_;
        }

        GAPP_API static GeneratorList* tls_generators() noexcept
        {
            static GeneratorList* tls_generators_ = new GeneratorList();
            return tls_generators_;
        }

        alignas(128) static thread_local RegisteredGenerator generator_;
    };


    /** The global pseudo-random number generator instance used in the algorithms. */
    inline constexpr ConcurrentXoroshiro128p prng;

} // namespace gapp::rng


/* IMPLEMENTATION */

namespace gapp::rng
{
    bool randomBool() noexcept
    {
        return bool(rng::prng() & 1u);
    }

    template<std::integral IntType>
    IntType randomInt(IntType lbound, IntType ubound) noexcept
    {
        GAPP_ASSERT(lbound <= ubound);

        return rng::uniform_int_distribution<IntType>{ lbound, ubound }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal() noexcept
    {
        return rng::generate_canonical<RealType>(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomReal(RealType lbound, RealType ubound) noexcept
    {
        GAPP_ASSERT(lbound <= ubound);

        return rng::uniform_real_distribution<RealType>{ lbound, ubound }(rng::prng);
    }

    template<std::floating_point RealType>
    RealType randomNormal(RealType mean, RealType std_dev) noexcept
    {
        GAPP_ASSERT(std_dev >= 0.0);

        return rng::normal_distribution<RealType>{ mean, std_dev }(rng::prng);
    }

    template<std::integral IntType>
    IntType randomBinomial(IntType n, double p) noexcept
    {
        GAPP_ASSERT(n >= 0);
        GAPP_ASSERT(0.0 <= p && p <= 1.0);

        return rng::binomial_distribution<IntType>{ n, p }(rng::prng);
    }

    template<std::integral IntType>
    void CachedRandomBinomial<IntType>::init(IntType n, double p) noexcept
    {
        dist_ = rng::binomial_distribution<IntType>(n, p);
    }

    template<std::integral IntType>
    IntType CachedRandomBinomial<IntType>::operator()(IntType n, double p) const noexcept
    {
        if (dist_.n() != n || dist_.p() != p)
        {
            return rng::binomial_distribution<IntType>{ n, p }(rng::prng);
        }

        return dist_(rng::prng);
    }

    template<std::ranges::random_access_range R>
    detail::size_type<R> randomIndex(const R& range) noexcept
    {
        GAPP_ASSERT(!range.empty());

        return rng::randomInt<detail::size_type<R>>(0, range.size() - 1);
    }

    template<std::forward_iterator Iter>
    Iter randomElement(Iter first, Iter last) noexcept
    {
        GAPP_ASSERT(std::distance(first, last) > 0);

        const auto max_offset = std::distance(first, last) - 1;
        const auto offset = rng::randomInt(0_pd, max_offset);

        return std::next(first, offset);
    }

    template<std::ranges::random_access_range Range, std::ranges::random_access_range Dist>
    std::ranges::range_reference_t<Range> randomElement(Range&& range, Dist&& cdf)
    {
        GAPP_ASSERT(!range.empty());
        GAPP_ASSERT(range.size() == cdf.size());

        return range[rng::sampleCdf(cdf)];
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
        const std::uint64_t range_len = detail::range_length(lbound, ubound);

        GAPP_ASSERT(ubound >= lbound);
        GAPP_ASSERT(range_len >= count);

        const bool select_many = (count >= std::uint64_t(0.6 * range_len));
        const bool huge_range  = (range_len >= 65536);

        if (huge_range) return rng::sampleUniqueSet(lbound, ubound, count);

        small_vector<IntType> numbers(count);

        thread_local detail::dynamic_bitset is_selected;
        is_selected.resize(range_len);
        is_selected.fill(select_many);

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

    template<std::ranges::random_access_range Range>
    size_t sampleCdf(const Range& cdf)
    {
        static_assert(std::is_floating_point_v<detail::value_t<Range>>);

        GAPP_ASSERT(!cdf.empty());
        GAPP_ASSERT(0.0 <= cdf.front());

        const auto threshold = rng::randomReal() * cdf.back(); // multiply by cdf.back() in case it's not exactly 1.0

        return std::distance(cdf.begin(), detail::lower_bound(cdf.begin(), cdf.end(), threshold));
    }

} // namespace gapp::rng

#endif // !GAPP_UTILITY_RNG_HPP
