/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATH_HPP
#define GA_UTILITY_MATH_HPP

#include "utility.hpp"
#include <vector>
#include <atomic>
#include <concepts>
#include <limits>
#include <cstdint>
#include <cstddef>

namespace genetic_algorithm::math
{
    class Tolerances
    {
    public:
        Tolerances() = delete;

        template<std::floating_point T = double>
        static T abs() noexcept { return T(absolute_tolerance.load(std::memory_order_acquire)); }

        template<std::floating_point T = double>
        static T eps() noexcept { return relative_tolerance_epsilons.load(std::memory_order_acquire) * std::numeric_limits<T>::epsilon(); }

    private:
        GA_API static std::atomic<double> absolute_tolerance;
        GA_API static std::atomic<unsigned> relative_tolerance_epsilons;

        friend class ScopedTolerances;
    };

    class ScopedTolerances
    {
    public:
        ScopedTolerances(unsigned num_epsilons, double abs) noexcept :
            old_abs_tol(Tolerances::absolute_tolerance.exchange(abs, std::memory_order_acq_rel)),
            old_eps_tol(Tolerances::relative_tolerance_epsilons.exchange(num_epsilons, std::memory_order_acq_rel))
        {}

        ~ScopedTolerances() noexcept
        {
            Tolerances::absolute_tolerance.store(old_abs_tol, std::memory_order_release);
            Tolerances::relative_tolerance_epsilons.store(old_eps_tol, std::memory_order_release);
        }

        ScopedTolerances(const ScopedTolerances&)            = delete;
        ScopedTolerances(ScopedTolerances&&)                 = delete;
        ScopedTolerances& operator=(const ScopedTolerances&) = delete;
        ScopedTolerances& operator=(ScopedTolerances&&)      = delete;

    private:
        double old_abs_tol;
        unsigned old_eps_tol;
    };
    

    template<typename T>
    inline constexpr T inf = std::numeric_limits<T>::infinity();

    using Point = std::vector<double>;

    using vector_iterator       = std::vector<double>::iterator;
    using const_vector_iterator = std::vector<double>::const_iterator;

    /* Comparison function for floating point numbers. Returns -1 if (lhs < rhs), +1 if (lhs > rhs), and 0 if (lhs == rhs). */
    template<std::floating_point T>
    constexpr std::int8_t floatCompare(T lhs, T rhs) noexcept;

    /* Equality comparison for floating point numbers. Returns true if lhs is approximately equal to rhs. */
    template<std::floating_point T>
    constexpr bool floatIsEqual(T lhs, T rhs) noexcept;

    /* Less than comparison for floating point numbers. Returns true if lhs is definitely less than rhs. */
    template<std::floating_point T>
    constexpr bool floatIsLess(T lhs, T rhs) noexcept;

    /* Less than comparison for fp numbers. Assumes that lhs is not greater than rhs. */
    template<std::floating_point T>
    constexpr bool floatIsLessAssumeNotGreater(T lhs, T rhs) noexcept;

    /* Greater than comparison for floating point numbers. Returns true if lhs is definitely greater than rhs. */
    template<std::floating_point T>
    constexpr bool floatIsGreater(T lhs, T rhs) noexcept;

    /* Less than or equal to comparison for floating point numbers. Returns true if lhs is less than or approximately equal to rhs. */
    template<std::floating_point T>
    constexpr bool floatIsLessEq(T lhs, T rhs) noexcept;

    /* Greater than or equal to comparison for floating point numbers. Returns true if lhs is greater than or approximately equal to rhs. */
    template<std::floating_point T>
    constexpr bool floatIsGreaterEq(T lhs, T rhs) noexcept;

    /* Equality comparison for fp vectors. Returns true if the elements of the vectors are approximately equal. */
    template<std::floating_point T>
    bool floatVecIsEqual(const std::vector<T>& lhs, const std::vector<T>& rhs) noexcept;

    /* Pareto comparison for fp vectors. Returns true if lhs is dominated by rhs (lhs < rhs), assuming maximization. */
    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs) noexcept;

    /* Pareto comparison for fp vectors. Returns -1 if (lhs < rhs), 1 if (lhs > rhs), and 0 if (lhs == rhs). */
    std::int8_t paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs) noexcept;
    
    /* Calculate the length of a vector. */
    double euclideanNorm(const std::vector<double>& vec) noexcept;

    /* Normalize the vector vec (divide by magnitude). */
    std::vector<double> normalizeVector(const std::vector<double>& vec);

    /* Normalize the vector vec (divide by magnitude). */
    std::vector<double> normalizeVector(std::vector<double>&& vec) noexcept;

    /* Calculate the square of the Euclidean distance between the vectors v1 and v2. */
    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2) noexcept;

    /* Calculate the square of the Euclidean distance between the vectors [first1, last1), [first2, first2 + last1 - first1). */
    double euclideanDistanceSq(const_vector_iterator first1, const_vector_iterator last1, const_vector_iterator first2) noexcept;

    /* Calculate the square of the perpendicular distance between a line (passing through the origin) and a point. */
    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point) noexcept;

    /* Calculate the square of the perpendicular distance between a line (passing through the origin) and a point. */
    double perpendicularDistanceSq(const_vector_iterator line_first, const_vector_iterator line_last, const_vector_iterator point_first) noexcept;

    /* Calculate the arithmetic mean of the values in vec. */
    double mean(const std::vector<double>& vec) noexcept;

    /* Calculate the standard deviation of the values in vec. */
    double stdDev(const std::vector<double>& vec) noexcept;

    /* Calculate the standard deviation of the values in vec. */
    double stdDev(const std::vector<double>& vec, double mean) noexcept;

    /* Computes the value of the function [ f(x) = integral sin(x)^n dx ] at x. */
    double integralSinPow(size_t exponent, double x) noexcept;

} // namespace genetic_algorithm::math


/* IMPLEMENTATION */

#include <algorithm>
#include <cmath>

namespace genetic_algorithm::math
{
    template<std::floating_point T>
    constexpr std::int8_t floatCompare(T lhs, T rhs) noexcept
    {
        // if lhs or rhs == NaN return unordered

        const T diff = lhs - rhs;
        const T scale = std::max(std::abs(lhs), std::abs(rhs));
        const T rel_tol = std::min(scale, std::numeric_limits<T>::max()) * Tolerances::eps<T>();
        const T tol = std::max(rel_tol, Tolerances::abs<T>());

        if (diff >  tol) return  1;  // lhs < rhs
        if (diff < -tol) return -1;  // lhs > rhs
        return 0;                    // lhs == rhs
    }

    template<std::floating_point T>
    constexpr bool floatIsEqual(T lhs, T rhs) noexcept
    {
        const T scale = std::max(std::abs(lhs), std::abs(rhs));

        if (scale == inf<T>) return lhs == rhs; // for infinities

        return std::abs(lhs - rhs) <= std::max(scale * Tolerances::eps<T>(), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsLess(T lhs, T rhs) noexcept
    {
        const T scale = std::max(std::abs(lhs), std::abs(rhs));

        if (scale == inf<T>) return lhs < rhs; // for infinities

        return (rhs - lhs) > std::max(scale * Tolerances::eps<T>(), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsLessAssumeNotGreater(T lhs, T rhs) noexcept
    {
        const T scale = std::abs(rhs);

        if (scale == inf<T>) return lhs < rhs; // for infinities

        return (rhs - lhs) > std::max(scale * Tolerances::eps<T>(), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsGreater(T lhs, T rhs) noexcept
    {
        const T scale = std::max(std::abs(lhs), std::abs(rhs));

        if (scale == inf<T>) return lhs > rhs; // for infinities

        return (lhs - rhs) > std::max(scale * Tolerances::eps<T>(), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsLessEq(T lhs, T rhs) noexcept
    {
        return !floatIsGreater(lhs, rhs);
    }

    template<std::floating_point T>
    constexpr bool floatIsGreaterEq(T lhs, T rhs) noexcept
    {
        return !floatIsLess(lhs, rhs);
    }

    template<std::floating_point T>
    bool floatVecIsEqual(const std::vector<T>& lhs, const std::vector<T>& rhs) noexcept
    {
        return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), floatIsEqual<T>);
    }

} // namespace genetic_algorithm::math

#endif // !GA_UTILITY_MATH_HPP