/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

#include "math.hpp"

#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    bool floatIsEqual(double lhs, double rhs, double eps)
    {
        assert(0.0 <= eps && eps <= 1.0);

        if (lhs == 0.0 || rhs == 0.0)
        {
            return std::abs(lhs - rhs) <= eps;
        }
        else
        {
            return std::abs(lhs - rhs) <= std::max(std::abs(lhs), std::abs(rhs)) * eps;
        }
    }

    bool floatIsLess(double lhs, double rhs, double eps)
    {
        assert(0.0 <= eps && eps <= 1.0);

        return (rhs - lhs) > std::max(std::abs(lhs), std::abs(rhs)) * eps;
    }

    bool floatIsEqual(const std::vector<double>& lhs, const std::vector<double>& rhs, double eps)
    {
        assert(0.0 <= eps && eps <= 1.0);
        assert(lhs.size() == rhs.size());

        return std::equal(lhs.begin(), lhs.end(), rhs.begin(), [eps](double lhs, double rhs) { return floatIsEqual(lhs, rhs, eps); });
    }

    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs, double eps)
    {
        assert(0.0 <= eps && eps <= 1.0);
        assert(lhs.size() == rhs.size());

        bool has_lower = false;
        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (floatIsLess(rhs[i], lhs[i], eps)) return false;
            if (floatIsLess(lhs[i], rhs[i], eps)) has_lower = true;
        }

        return has_lower;
    }

    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2)
    {
        assert(v1.size() == v2.size());

        return std::transform_reduce(v1.begin(), v1.end(), v2.begin(), 0.0,
                                     std::plus<double>(),
                                     [](double lhs, double rhs) { return std::pow(lhs - rhs, 2); });
    }

    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point)
    {
        assert(line.size() == point.size());
        assert(!line.empty());

        double k = std::inner_product(line.begin(), line.end(), point.begin(), 0.0) /
                   std::inner_product(line.begin(), line.end(), line.begin(), 0.0);

        return std::transform_reduce(point.begin(), point.end(), line.begin(), 0.0,
                                     std::plus<double>(),
                                     [k](double p, double l) { return std::pow(p - k * l, 2); });
    }
}