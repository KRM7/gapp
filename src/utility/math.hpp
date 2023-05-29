/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATH_HPP
#define GA_UTILITY_MATH_HPP

#include "utility.hpp"
#include <vector>
#include <span>
#include <atomic>
#include <concepts>
#include <limits>
#include <cstdint>
#include <cstddef>

/** Math utility classes and functions. */
namespace gapp::math
{
    /**
    * This class contains the global absolute and relative tolerance values used
    * for comparing floating-point values in the GAs.
    * 
    * New tolerances can be set using the ScopedTolerances class.
    */
    class Tolerances
    {
    public:
        Tolerances() = delete;

        /** @returns The current absolute tolerance used for floating-point comparisons. */
        template<std::floating_point T = double>
        static T abs() noexcept { return T(absolute_tolerance.load(std::memory_order_acquire)); }

        /** @returns The current relative tolerance used for floating-point comparisons. */
        template<std::floating_point T = double>
        static T eps() noexcept { return relative_tolerance_epsilons.load(std::memory_order_acquire) * std::numeric_limits<T>::epsilon(); }

    private:
        GA_API static std::atomic<double> absolute_tolerance;
        GA_API static std::atomic<unsigned> relative_tolerance_epsilons;

        friend class ScopedTolerances;
    };

    /**
    * This class can be used to set the absolute and relative tolerances used for
    * comparing floating-point values in the GAs.
    * 
    * When an instance of the class is created, the tolerance values are set to the
    * values specified by the parameters of the constructor, and these new tolerance
    * values will be used for floating-point comparisons until the instance of the class is destroyed.
    * The tolerances are reset to their old values when the instance is destroyed.
    */
    class [[nodiscard]] ScopedTolerances
    {
    public:
        /**
        * Create an instance of the class, setting new values for the tolerances
        * used for floating-point comparisons.
        * 
        * @param num_epsilons The number of epsilons to use as the relative tolerance in the comparisons.
        * @param abs The absolute tolerance value used for the comparisons.
        */
        ScopedTolerances(unsigned num_epsilons, double abs) noexcept :
            old_abs_tol(Tolerances::absolute_tolerance.exchange(abs, std::memory_order_acq_rel)),
            old_eps_tol(Tolerances::relative_tolerance_epsilons.exchange(num_epsilons, std::memory_order_acq_rel))
        {}

        /** Reset the tolerances to their previous values. */
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

    template<typename T>
    inline constexpr T small = std::numeric_limits<T>::min();

    template<typename T>
    inline constexpr T large = std::numeric_limits<T>::max();


    using Point = std::vector<double>;


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

    /* Equality comparison for fp vectors. Returns true if the elements of the ranges are approximately equal. */
    template<std::floating_point T>
    bool floatVecIsEqual(std::span<const T> lhs, std::span<const T> rhs) noexcept;

    /* Pareto comparison for fp ranges. Returns true if lhs is dominated by rhs (lhs < rhs), assuming maximization. */
    bool paretoCompareLess(std::span<const double> lhs, std::span<const double> rhs) noexcept;

    /* Pareto comparison for fp ranges. Returns -1 if (lhs < rhs), 1 if (lhs > rhs), and 0 if (lhs == rhs). */
    std::int8_t paretoCompare(std::span<const double> lhs, std::span<const double> rhs) noexcept;
    
    /* Calculate the length of a vector. */
    double euclideanNorm(std::span<const double> vec) noexcept;

    /* Normalize the vector vec (divide by magnitude). */
    std::vector<double> normalizeVector(const std::vector<double>& vec);

    /* Normalize the vector vec (divide by magnitude). */
    std::vector<double> normalizeVector(std::vector<double>&& vec) noexcept;

    /* Calculate the square of the Euclidean distance between the vectors v1 and v2. */
    double euclideanDistanceSq(std::span<const double> v1, std::span<const double> v2) noexcept;

    /* Calculate the square of the perpendicular distance between a line (passing through the origin) and a point. */
    double perpendicularDistanceSq(std::span<const double> line, std::span<const double> point) noexcept;

    /* Calculate the volume of a hyperrectangle given by 2 points. */
    double volumeBetween(std::span<const double> p1, std::span<const double> p2) noexcept;


    /* Calculate the arithmetic mean of the values in vec. */
    double mean(std::span<const double> vec) noexcept;

    /* Calculate the standard deviation of the values in vec. */
    double stdDev(std::span<const double> vec) noexcept;

    /* Calculate the standard deviation of the values in vec. */
    double stdDev(std::span<const double> vec, double mean) noexcept;


    /* Computes the value of the function [ f(x) = integral sin(x)^n dx ] at x. */
    double integralSinPow(size_t exponent, double x) noexcept;

} // namespace gapp::math


/* IMPLEMENTATION */

#include <algorithm>
#include <cmath>

namespace gapp::math
{
    template<std::floating_point T>
    constexpr std::int8_t floatCompare(T lhs, T rhs) noexcept
    {
        GA_ASSERT(!std::isnan(lhs) && !std::isnan(rhs));

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
    bool floatVecIsEqual(std::span<const T> lhs, std::span<const T> rhs) noexcept
    {
        return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), floatIsEqual<T>);
    }

} // namespace gapp::math

#endif // !GA_UTILITY_MATH_HPP