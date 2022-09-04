/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "math.hpp"
#include "functional.hpp"
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::math
{
    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs) noexcept
    {
        return paretoCompareLess(lhs, rhs, 0);
    }

    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs, size_t first) noexcept
    {
        assert(lhs.size() == rhs.size());

        for (size_t i = first; i < lhs.size(); i++)
        {
            if (floatIsLess(rhs[i], lhs[i])) return false;
        }
        for (size_t i = first; i < lhs.size(); i++)
        {
            if (floatIsLessAssumeNotGreater(lhs[i], rhs[i])) return true;
        }

        return false;
    }

    std::int8_t paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs) noexcept
    {
        return paretoCompare(lhs, rhs, 0);
    }

    std::int8_t paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs, size_t first) noexcept
    {
        assert(lhs.size() == rhs.size());

        std::int8_t lhs_has_lower = 0;
        std::int8_t rhs_has_lower = 0;
        for (size_t i = first; i < lhs.size(); i++)
        {
            auto comp = floatCompare(lhs[i], rhs[i]);
            if (comp < 0)
            {
                if (rhs_has_lower) return 0;
                lhs_has_lower = 1;
            }
            else if (comp > 0)
            {
                if (lhs_has_lower) return 0;
                rhs_has_lower = 1;
            }
        }

        return rhs_has_lower - lhs_has_lower;
    }

    double euclideanNorm(const std::vector<double>& vec) noexcept
    {
        return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0));
    }

    std::vector<double> normalizeVector(const std::vector<double>& vec) noexcept
    {
        const double mag = euclideanNorm(vec);

        return detail::map(vec, detail::divide_by(mag));
    }

    std::vector<double> normalizeVector(std::vector<double>&& vec) noexcept
    {
        const double mag = euclideanNorm(vec);

        std::transform(vec.begin(), vec.end(), vec.begin(), detail::divide_by(mag));

        return std::move(vec);
    }

    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2) noexcept
    {
        assert(v1.size() == v2.size());

        return std::transform_reduce(v1.begin(), v1.end(), v2.begin(), 0.0,
                                     std::plus{},
                                     [](double lhs, double rhs) noexcept { return std::pow(lhs - rhs, 2); });
    }

    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point) noexcept
    {
        assert(line.size() == point.size());
        assert(!line.empty());

        double k = std::inner_product(line.begin(), line.end(), point.begin(), 0.0) /
                   std::inner_product(line.begin(), line.end(), line.begin(),  0.0);

        return std::transform_reduce(point.begin(), point.end(), line.begin(), 0.0,
                                     std::plus{},
                                     [k](double p, double l) noexcept { return std::pow(p - k * l, 2); });
    }

    double mean(const std::vector<double>& vec) noexcept
    {
        assert(!vec.empty());

        return std::transform_reduce(vec.begin(), vec.end(), 0.0, std::plus{}, detail::divide_by(double(vec.size())));
    }

    double stdDev(const std::vector<double>& vec) noexcept
    {
        return stdDev(vec, mean(vec));
    }

    double stdDev(const std::vector<double>& vec, double mean) noexcept
    {
        assert(!vec.empty());

        if (vec.size() == 1) return 0.0;

        auto var = std::transform_reduce(vec.begin(), vec.end(), 0.0, std::plus{},
        [mean, n = 1.0 / (vec.size() - 1)](double val) noexcept
        {
            return std::pow(val - mean, 2) * n;
        });

        return std::sqrt(var);
    }

} // namespace genetic_algorithm::math