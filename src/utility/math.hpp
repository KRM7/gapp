/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATH_HPP
#define GA_UTILITY_MATH_HPP

#include "utility.hpp"
#include <vector>
#include <concepts>
#include <cstdint>

namespace genetic_algorithm::math
{
    class Tolerances
    {
    public:
        Tolerances() = delete;

        template<std::floating_point T = double>
        static T abs() noexcept { return T(absolute_tolerance); }

        template<std::floating_point T = double>
        static T eps() noexcept { return relative_tolerance_epsilons * std::numeric_limits<T>::epsilon(); }

        static void abs(double tolerance) { assert(tolerance >= 0.0); absolute_tolerance = tolerance; }
        static void eps(unsigned n) noexcept { relative_tolerance_epsilons = n; }

    private:
        static double absolute_tolerance;
        static unsigned relative_tolerance_epsilons;

        friend class LocalTolerances;
    };

    class LocalTolerances
    {
    public:
        LocalTolerances(unsigned num_epsilons, double abs) :
            old_abs_tol(Tolerances::absolute_tolerance),
            old_eps_tol(Tolerances::relative_tolerance_epsilons)
        {
            Tolerances::abs(abs);
            Tolerances::eps(num_epsilons);
        }

        ~LocalTolerances()
        {
            Tolerances::abs(old_abs_tol);
            Tolerances::eps(old_eps_tol);
        }

        LocalTolerances(const LocalTolerances&)            = delete;
        LocalTolerances(LocalTolerances&&)                 = delete;
        LocalTolerances& operator=(const LocalTolerances&) = delete;
        LocalTolerances& operator=(LocalTolerances&&)      = delete;

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

    /*  */
    double integralSinPow(size_t exponent, double x) noexcept;

} // namespace genetic_algorithm::math


/* IMPLEMENTATION */

#include <algorithm>
#include <cmath>
#include <cassert>

namespace genetic_algorithm::math
{
    template<std::floating_point T>
    constexpr std::int8_t floatCompare(T lhs, T rhs) noexcept
    {
        if (lhs == rhs) return 0; // for infinities

        const T diff = lhs - rhs;
        const T rel_tol = std::max(std::abs(lhs), std::abs(rhs)) * Tolerances::eps<T>();
        const T tol = std::max(rel_tol, Tolerances::abs<T>());

        if (diff >= tol)  return  1;  // lhs < rhs
        if (diff <= -tol) return -1;  // lhs > rhs
        return 0;                     // lhs == rhs
    }

    template<std::floating_point T>
    constexpr bool floatIsEqual(T lhs, T rhs) noexcept
    {
        if (lhs == rhs) return true; // for infinities

        const T rel_tol = std::max(std::abs(lhs), std::abs(rhs)) * Tolerances::eps<T>();

        return std::abs(lhs - rhs) < std::max(rel_tol, Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsLess(T lhs, T rhs) noexcept
    {
        const T rel_tol = std::max(std::abs(lhs), std::abs(rhs)) * Tolerances::eps<T>();

        return (rhs - lhs) >= std::max(rel_tol, Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsLessAssumeNotGreater(T lhs, T rhs) noexcept
    {
        const T rel_tol = std::abs(rhs) * Tolerances::eps<T>();

        return (rhs - lhs) >= std::max(rel_tol, Tolerances::abs<T>());
    }

    template<std::floating_point T>
    constexpr bool floatIsGreater(T lhs, T rhs) noexcept
    {
        const T rel_tol = std::max(std::abs(lhs), std::abs(rhs)) * Tolerances::eps<T>();

        return (lhs - rhs) >= std::max(rel_tol, Tolerances::abs<T>());
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