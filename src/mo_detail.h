/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
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

/**
* This file contains some utility functions for the NSGA-III algorithm.
*/

#ifndef GA_MO_DETAIL_H
#define GA_MO_DETAIL_H

#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::detail
{
    /* Return true if lhs is dominated by rhs (lhs < rhs) assuming maximization. */
    inline bool paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs);

    /* Calculate the square of the Euclidean distance between the vectors v1 and v2. */
    inline double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2);

    /* Calculate the square of the perpendicular distance between the line ref and the point p. */
    inline double perpendicularDistanceSq(const std::vector<double>& ref, const std::vector<double>& p);

    /* Find the index and distance of the closest reference line to the point p. */
    inline std::pair<size_t, double> findClosestRef(const std::vector<std::vector<double>>& refs, const std::vector<double>& p);

    /* Achievement scalarization function. */
    inline double ASF(const std::vector<double>& f, const std::vector<double>& z, const std::vector<double>& w);

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <algorithm>
#include <utility>
#include <cmath>
#include <cassert>

namespace genetic_algorithm::detail
{
    bool paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs)
    {
        assert(lhs.size() == rhs.size());

        bool has_lower = false;
        for (size_t i = 0; i < lhs.size(); i++)
        {
            if (lhs[i] > rhs[i]) return false;
            if (lhs[i] < rhs[i]) has_lower = true;
        }
        return has_lower;
    }

    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2)
    {
        assert(v1.size() == v2.size());

        double d = 0.0;
        for (size_t i = 0; i < v1.size(); i++)
        {
            d += (v1[i] - v2[i]) * (v1[i] - v2[i]);
        }

        return d;
    }

    double perpendicularDistanceSq(const std::vector<double>& ref, const std::vector<double>& p)
    {
        assert(ref.size() == p.size());

        double num = 0.0, den = 0.0;
        for (size_t i = 0; i < ref.size(); i++)
        {
            num += ref[i] * p[i];
            den += ref[i] * ref[i];
        }
        double k = num / den;

        double dist = 0.0;
        for (size_t i = 0; i < ref.size(); i++)
        {
            dist += (p[i] - k * ref[i]) * (p[i] - k * ref[i]);
        }

        return dist;
    }

    std::pair<size_t, double> findClosestRef(const std::vector<std::vector<double>>& refs, const std::vector<double>& p)
    {
        size_t argmin = 0;
        double dmin = perpendicularDistanceSq(refs[0], p);
        for (size_t i = 1; i < refs.size(); i++)
        {
            double d = perpendicularDistanceSq(refs[i], p);
            if (d < dmin)
            {
                dmin = d;
                argmin = i;
            }
        }

        return std::make_pair(argmin, dmin);
    }

    double ASF(const std::vector<double>& f, const std::vector<double>& z, const std::vector<double>& w)
    {
        assert(!f.empty());
        assert(f.size() == z.size() && f.size() == w.size());

        double dmax = std::abs(f[0] - z[0]) / w[0];
        for (size_t j = 1; j < f.size(); j++)
        {
            dmax = std::max(dmax, std::abs(f[j] - z[j]) / w[j]);
        }

        return dmax;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_MO_DETAIL_H