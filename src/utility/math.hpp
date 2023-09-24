/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATH_HPP
#define GA_UTILITY_MATH_HPP

#include "bounded_value.hpp"
#include "utility.hpp"
#include <vector>
#include <span>
#include <compare>
#include <concepts>
#include <limits>
#include <cstdint>
#include <cstddef>

/** Math utility classes and functions. */
namespace gapp::math
{
    template<typename T>
    inline constexpr T inf = std::numeric_limits<T>::infinity();

    template<typename T>
    inline constexpr T eps = std::numeric_limits<T>::epsilon();

    template<typename T>
    inline constexpr T small = std::numeric_limits<T>::min();

    template<typename T>
    inline constexpr T large = std::numeric_limits<T>::max();

    using Point = std::vector<double>;


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
        static T abs() noexcept { return T(absolute_tolerance); }

        /** @returns The current relative tolerance used for floating-point comparisons around @p at. */
        template<std::floating_point T = double>
        static T rel(T at) noexcept { return relative_tolerance * at; }

    private:
        GAPP_API inline constinit static double absolute_tolerance = 1E-12;
        GAPP_API inline constinit static double relative_tolerance = 10 * eps<double>;

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
    * 
    * @warning
    *   Creating an instance of this class modifies the global floating point tolerance
    *   values. This means that this class shouldn't be instantiated on multiple threads
    *   at once.
    */
    class [[nodiscard]] ScopedTolerances
    {
    public:
        /**
        * Create an instance of the class, setting new values for the tolerances
        * used for floating-point comparisons.
        * 
        * @param abs The absolute tolerance value that will be used for the comparisons. Can't be negative.
        * @param rel The relative tolerance value around 1.0 that will be used for the comparisons. Can't be negative.
        */
        ScopedTolerances(NonNegative<double> abs, NonNegative<double> rel) noexcept :
            old_absolute_tolerance(std::exchange(Tolerances::absolute_tolerance, abs)),
            old_relative_tolerance(std::exchange(Tolerances::relative_tolerance, rel))
        {}

        /** Reset the tolerances to their previous values. */
        ~ScopedTolerances() noexcept
        {
            Tolerances::absolute_tolerance = old_absolute_tolerance;
            Tolerances::relative_tolerance = old_relative_tolerance;
        }

        ScopedTolerances(const ScopedTolerances&)            = delete;
        ScopedTolerances(ScopedTolerances&&)                 = delete;
        ScopedTolerances& operator=(const ScopedTolerances&) = delete;
        ScopedTolerances& operator=(ScopedTolerances&&)      = delete;

    private:
        double old_absolute_tolerance;
        double old_relative_tolerance;
    };


    /* Three-way comparison function for floating point numbers. */
    template<std::floating_point T>
    std::weak_ordering floatCompare(T lhs, T rhs) noexcept;

    /* Equality comparison for floating point numbers. Returns true if lhs is approximately equal to rhs. */
    template<std::floating_point T>
    bool floatIsEqual(T lhs, T rhs) noexcept;

    /* Less than comparison for floating point numbers. Returns true if lhs is definitely less than rhs. */
    template<std::floating_point T>
    bool floatIsLess(T lhs, T rhs) noexcept;

    /* Less than comparison for fp numbers. Assumes that lhs is not greater than rhs. */
    template<std::floating_point T>
    bool floatIsLessAssumeNotGreater(T lhs, T rhs) noexcept;

    /* Greater than comparison for floating point numbers. Returns true if lhs is definitely greater than rhs. */
    template<std::floating_point T>
    bool floatIsGreater(T lhs, T rhs) noexcept;

    /* Less than or equal to comparison for floating point numbers. Returns true if lhs is less than or approximately equal to rhs. */
    template<std::floating_point T>
    bool floatIsLessEq(T lhs, T rhs) noexcept;

    /* Greater than or equal to comparison for floating point numbers. Returns true if lhs is greater than or approximately equal to rhs. */
    template<std::floating_point T>
    bool floatIsGreaterEq(T lhs, T rhs) noexcept;

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
    std::weak_ordering floatCompare(T lhs, T rhs) noexcept
    {
        GAPP_ASSERT(!std::isnan(lhs) && !std::isnan(rhs));

        const T diff = lhs - rhs;
        const T scale = std::min(std::max(std::abs(lhs), std::abs(rhs)), std::numeric_limits<T>::max());
        const T tol = std::max(Tolerances::rel<T>(scale), Tolerances::abs<T>());

        if (diff >  tol) return std::weak_ordering::greater;
        if (diff < -tol) return std::weak_ordering::less;
        return std::weak_ordering::equivalent;
    }

    template<std::floating_point T>
    bool floatIsEqual(T lhs, T rhs) noexcept
    {
        const T scale = std::max(std::abs(lhs), std::abs(rhs));

        if (scale == inf<T>) return lhs == rhs; // for infinities

        return std::abs(lhs - rhs) <= std::max(Tolerances::rel<T>(scale), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    bool floatIsLess(T lhs, T rhs) noexcept
    {
        const T scale = std::max(std::abs(lhs), std::abs(rhs));

        if (scale == inf<T>) return lhs < rhs; // for infinities

        return (rhs - lhs) > std::max(Tolerances::rel<T>(scale), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    bool floatIsLessAssumeNotGreater(T lhs, T rhs) noexcept
    {
        const T scale = std::abs(rhs);

        if (scale == inf<T>) return lhs < rhs; // for infinities

        return (rhs - lhs) > std::max(Tolerances::rel<T>(scale), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    bool floatIsGreater(T lhs, T rhs) noexcept
    {
        const T scale = std::max(std::abs(lhs), std::abs(rhs));

        if (scale == inf<T>) return lhs > rhs; // for infinities

        return (lhs - rhs) > std::max(Tolerances::rel<T>(scale), Tolerances::abs<T>());
    }

    template<std::floating_point T>
    bool floatIsLessEq(T lhs, T rhs) noexcept
    {
        return !floatIsGreater(lhs, rhs);
    }

    template<std::floating_point T>
    bool floatIsGreaterEq(T lhs, T rhs) noexcept
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