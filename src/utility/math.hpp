/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MATH_HPP
#define GA_MATH_HPP

#include "utility.hpp"
#include <vector>
#include <concepts>

namespace genetic_algorithm::detail
{
    /* Equality comparison for floating point numbers. Returns true if lhs is approximately equal to rhs. */
    template<std::floating_point T>
    bool floatIsEqual(T lhs, T rhs) noexcept;

    /* Less than comparison for floating point numbers. Returns true if lhs is definitely less than rhs. */
    template<std::floating_point T>
    bool floatIsLess(T lhs, T rhs) noexcept;

    /* Greater than comparison for floating point numbers. Returns true if lhs is definitely greater than rhs. */
    template<std::floating_point T>
    bool floatIsGreater(T lhs, T rhs) noexcept;

    /* Less than or equal to comparison for floating point numbers. Returns true if lhs is less than or approximately equal to rhs. */
    template<std::floating_point T>
    bool floatIsLessEq(T lhs, T rhs) noexcept;

    /* Greater than or equal to comparison for floating point numbers. Returns true if lhs is greater than or approximately equal to rhs. */
    template<std::floating_point T>
    bool floatIsGreaterEq(T lhs, T rhs) noexcept;

    /* Equality comparison for fp vectors. Returns true if the elements of the vectors are approximately equal. */
    template<std::floating_point T>
    bool floatVecIsEqual(const std::vector<T>& lhs, const std::vector<T>& rhs) noexcept;

    /* Pareto comparison for fp vectors. Returns true if lhs is dominated by rhs (lhs < rhs), assuming maximization. */
    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs) noexcept;

    /* Pareto comparison starting at idx = first. */
    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs, size_t first) noexcept;

    /* Calculate the square of the Euclidean distance between the vectors v1 and v2. */
    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2) noexcept;

    /* Calculate the square of the perpendicular distance between a line and a point. */
    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point) noexcept;

    /* Calculate the arithmetic mean of the values in vec. */
    double mean(const std::vector<double>& vec) noexcept;

    /* Calculate the standard deviation of the values in vec. */
    double stdDev(const std::vector<double>& vec) noexcept;

    /* Calculate the standard deviation of the values in vec. */
    double stdDev(const std::vector<double>& vec, double mean) noexcept;

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <algorithm>
#include <cmath>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<std::floating_point T>
    bool floatIsEqual(T lhs, T rhs) noexcept
    {
        return std::abs(lhs - rhs) <= std::max({ std::abs(lhs), std::abs(rhs) }) * GA_EPSILON;
    }

    template<std::floating_point T>
    bool floatIsLess(T lhs, T rhs) noexcept
    {
        return (rhs - lhs) > std::max({ std::abs(lhs), std::abs(rhs) }) * GA_EPSILON;
    }

    template<std::floating_point T>
    bool floatIsGreater(T lhs, T rhs) noexcept
    {
        return (lhs - rhs) > std::max({ std::abs(lhs), std::abs(rhs) }) * GA_EPSILON;
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
    bool floatVecIsEqual(const std::vector<T>& lhs, const std::vector<T>& rhs) noexcept
    {
        return (lhs.size() == rhs.size()) &&
                std::equal(lhs.begin(), lhs.end(), rhs.begin(), floatIsEqual<T>);
    }

} // namespace genetic_algorithm::detail

#endif // !GA_MATH_HPP