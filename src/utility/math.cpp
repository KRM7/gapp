/* Copyright (c) 2022 Kriszti�n Rug�si. Subject to the MIT License. */

#include "math.hpp"
#include "algorithm.hpp"
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    constexpr double square(double n) noexcept
    {
        return n * n;
    }

    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs) noexcept
    {
        return paretoCompareLess(lhs, rhs, 0);
    }

    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs, size_t first) noexcept
    {
        assert(lhs.size() == rhs.size());

        bool has_lower = false;
        for (size_t i = first; i < lhs.size(); i++)
        {
            if (floatIsLess(rhs[i], lhs[i])) return false;
            if (floatIsLess(lhs[i], rhs[i])) has_lower = true;
        }

        return has_lower;
    }

        for (size_t i = first; i < lhs.size(); i++)
        {
        }

    }

    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2) noexcept
    {
        assert(v1.size() == v2.size());

        return std::transform_reduce(v1.begin(), v1.end(), v2.begin(), 0.0,
                                     std::plus<double>{},
                                     compose(std::minus<double>{}, square));
    }

    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point) noexcept
    {
        assert(line.size() == point.size());
        assert(!line.empty());

        double k = std::inner_product(line.begin(), line.end(), point.begin(), 0.0) /
                   std::inner_product(line.begin(), line.end(), line.begin(), 0.0);

        return std::transform_reduce(point.begin(), point.end(), line.begin(), 0.0,
                                     std::plus<double>{},
                                     [k](double p, double l) noexcept { return square(p - k * l); });
    }

    double mean(const std::vector<double>& vec) noexcept
    {
        assert(!vec.empty());

        return std::accumulate(vec.begin(), vec.end(), 0.0,
        [n = vec.size()](double acc, double val) noexcept
        {
            return acc + val / n;
        });
    }

    double stdDev(const std::vector<double>& vec) noexcept
    {
        return stdDev(vec, mean(vec));
    }

    double stdDev(const std::vector<double>& vec, double mean) noexcept
    {
        assert(!vec.empty());

        if (vec.size() == 1) return 0.0;

        auto var = std::accumulate(vec.begin(), vec.end(), 0.0L,
        [mean, n = vec.size() - 1](long double acc, double val) noexcept
        {
            return acc + square(val - mean) / n;
        });

        return double(std::sqrt(var));
    }
}