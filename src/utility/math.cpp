/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "math.hpp"
#include "functional.hpp"
#include "utility.hpp"
#include <algorithm>
#include <numeric>
#include <span>
#include <functional>
#include <iterator>
#include <type_traits>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace gapp::math
{
    bool paretoCompareLess(std::span<const double> lhs, std::span<const double> rhs) noexcept
    {
        GAPP_ASSERT(lhs.size() == rhs.size());

        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (floatIsLess(rhs[i], lhs[i])) return false;
        }
        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (floatIsLessAssumeNotGreater(lhs[i], rhs[i])) return true;
        }

        return false;
    }

    std::int8_t paretoCompare(std::span<const double> lhs, std::span<const double> rhs) noexcept
    {
        GAPP_ASSERT(lhs.size() == rhs.size());

        std::int8_t lhs_has_lower = 0;
        std::int8_t rhs_has_lower = 0;

        for (size_t i = 0; i < lhs.size(); i++)
        {
            const auto comp = floatCompare(lhs[i], rhs[i]);
            if (comp < 0)       // NOLINT(*use-nullptr)
            {
                if (rhs_has_lower) return 0;
                lhs_has_lower = 1;
            }
            else if (comp > 0)  // NOLINT(*use-nullptr)
            {
                if (lhs_has_lower) return 0;
                rhs_has_lower = 1;
            }
        }

        return std::int8_t(rhs_has_lower - lhs_has_lower);
    }

    double euclideanNorm(std::span<const double> vec) noexcept
    {
        return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0));
    }

    void normalizeVector(std::span<double> vec) noexcept
    {
        const double norm = euclideanNorm(vec);

        std::transform(vec.begin(), vec.end(), vec.begin(), detail::divide_by(norm));
    }

    double euclideanDistanceSq(std::span<const double> v1, std::span<const double> v2) noexcept
    {
        GAPP_ASSERT(v1.size() == v2.size());

        return std::transform_reduce(v1.begin(), v1.end(), v2.begin(), 0.0,
                                     std::plus{},
                                     [](double lhs, double rhs) noexcept { return std::pow(lhs - rhs, 2); });
    }

    double perpendicularDistanceSq(std::span<const double> line, std::span<const double> point) noexcept
    {
        GAPP_ASSERT(line.size() == point.size());

        double k = std::inner_product(line.begin(), line.end(), point.begin(), 0.0) /
                   std::inner_product(line.begin(), line.end(), line.begin(), 0.0);

        return std::transform_reduce(line.begin(), line.end(), point.begin(), 0.0,
                                     std::plus{},
                                     [k](double l, double p) noexcept { return std::pow(p - k * l, 2); });
    }

    double volumeBetween(std::span<const double> p1, std::span<const double> p2) noexcept
    {
        GAPP_ASSERT(p1.size() == p2.size());

        return std::abs(std::transform_reduce(p1.begin(), p1.end(), p2.begin(), 1.0, std::multiplies{}, std::minus{}));
    }

    double mean(std::span<const double> vec) noexcept
    {
        GAPP_ASSERT(!vec.empty());

        return std::transform_reduce(vec.begin(), vec.end(), 0.0, std::plus{}, detail::divide_by(vec.size()));
    }

    double stdDev(std::span<const double> vec) noexcept
    {
        return stdDev(vec, mean(vec));
    }

    double stdDev(std::span<const double> vec, double mean) noexcept
    {
        GAPP_ASSERT(!vec.empty());

        if (vec.size() == 1) return 0.0;

        const double var = std::transform_reduce(vec.begin(), vec.end(), 0.0, std::plus{},
        [mean, n = 1.0 / std::sqrt(vec.size())](double val) noexcept
        {
            return std::pow(n * (val - mean), 2);
        });

        return std::sqrt(var);
    }

    double integralSinPow(size_t exponent, double x) noexcept
    {
        double integral = (exponent % 2) ? -std::cos(x) : x;

        for (size_t n = 2 + exponent % 2; n <= exponent; n += 2)
        {
            const double mult = (n - 1.0) / n;
            const double plus = -std::cos(x) * std::pow(std::sin(x), n - 1) / n;

            integral = mult * integral + plus;
        }

        return integral;
    }

} // namespace gapp::math
