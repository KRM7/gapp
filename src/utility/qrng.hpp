/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_QRNG_HPP
#define GA_UTILITY_QRNG_HPP

#include "bounded_value.hpp"
#include <vector>
#include <concepts>
#include <cstddef>

namespace gapp::rng
{
    /**
    * Quasi-random number generator based on: http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/.
    * Generates points in the unit hypercube in dim dimensions.
    */
    template<std::floating_point RealType = double>
    class QuasiRandom
    {
    public:
        using result_type = std::vector<RealType>;
        using state_type  = std::vector<RealType>;
        using size_type   = std::size_t;

        /** Create a quasi-random number generator in @p dim dimensions. */
        explicit QuasiRandom(size_type dim, NonNegative<RealType> seed = 0.5);

        /** Generate the next quasi-random point of the sequence. */
        result_type operator()();

        /** Discard the next @p n points of the sequence. */
        void discard(size_type n = 1);

        /** Reset the generator's state using a specified seed. */
        void reset(NonNegative<RealType> new_seed = 0.5);

        /** @returns The generator's number of dimensions. */
        [[nodiscard]]
        size_type dim() const noexcept;

    private:
        size_type dim_;         /* The dimension of the generated points of the sequence. */
        RealType seed_;         /* The seed used. */
        state_type alpha_;      /* The initial point of the sequence (without considering the seed). */
        state_type point_;      /* The current/last point generated in the sequence. */

        /* Approximate the generalized golden ratio in dim dimensions. */
        static constexpr RealType phi(size_t dim, size_t n = 30) noexcept;
    };

} // namespace gapp::rng


/* IMPLEMENTATION */

#include <cmath>
#include <algorithm>
#include "utility.hpp"

namespace gapp::rng
{
    template<std::floating_point T>
    QuasiRandom<T>::QuasiRandom(size_type dim, NonNegative<T> seed) :
        dim_(dim), seed_(seed), alpha_(dim, 0.0), point_(dim, 0.0)
    {
        const T phid = phi(dim);
        T phidpow = phid;

        for (size_type i = 0; i < dim; i++)
        {
            alpha_[i] = 1.0 / phidpow;
            point_[i] = seed_;
            phidpow *= phid;
        }
    }

    template<std::floating_point T>
    inline auto QuasiRandom<T>::operator()() -> result_type
    {
        for (size_type i = 0; i < point_.size(); i++)
        {
            point_[i] += alpha_[i];
            point_[i] -= size_type(point_[i]);
        }

        return point_;
    }

    template<std::floating_point RealType>
    inline void QuasiRandom<RealType>::discard(size_type n)
    {
        while (n--) (void)operator()();
    }

    template<std::floating_point RealType>
    inline void QuasiRandom<RealType>::reset(NonNegative<RealType> new_seed)
    {
        seed_ = new_seed;
        std::fill(point_.begin(), point_.end(), seed_);
    }

    template<std::floating_point RealType>
    inline auto QuasiRandom<RealType>::dim() const noexcept -> size_type
    {
        return dim_;
    }

    template<std::floating_point RealType>
    constexpr inline RealType QuasiRandom<RealType>::phi(size_type dim, size_t n) noexcept
    {
        RealType phid = 1.0;
        while (n--)
        {
            phid = std::pow(1.0 + phid, 1.0 / (dim + 1.0));
        }

        return phid;
    }

} // namespace gapp::rng

#endif // !GA_UTILITY_QRNG_HPP