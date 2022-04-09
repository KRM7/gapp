/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

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

    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs)
    {
        assert(0.0 <= eps && eps <= 1.0);
        assert(lhs.size() == rhs.size());

        bool has_lower = false;
        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (floatIsLess(rhs[i], lhs[i])) return false;
            if (floatIsLess(lhs[i], rhs[i])) has_lower = true;
        }

        return has_lower;
    }

    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2)
    {
        assert(v1.size() == v2.size());

        return std::transform_reduce(v1.begin(), v1.end(), v2.begin(), 0.0,
                                     std::plus<double>{},
                                     compose(std::minus<double>{}, square));
    }

    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point)
    {
        assert(line.size() == point.size());
        assert(!line.empty());

        double k = std::inner_product(line.begin(), line.end(), point.begin(), 0.0) /
                   std::inner_product(line.begin(), line.end(), line.begin(), 0.0);

        return std::transform_reduce(point.begin(), point.end(), line.begin(), 0.0,
                                     std::plus<double>(),
                                     [k](double p, double l) { return square(p - k * l); });
    }

    double mean(const std::vector<double>& vec) noexcept
    {
        assert(!vec.empty());

        return std::accumulate(vec.begin(), vec.end(), 0.0,
        [n = vec.size()](double acc, double val)
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
        [mean, n = vec.size()](long double acc, double val)
        {
            return acc + square(val - mean) / (n - 1.0);
        });

        return double(std::sqrt(var));
    }
}