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

/*
* This file contains the functions used for generating the reference points
* in the NSGA-III algorithm.
*/

#ifndef GA_REFERENCE_POINTS_H
#define GA_REFERENCE_POINTS_H

#include <vector>
#include <cstddef>

namespace genetic_algorithm::detail
{
    /* Sample a point from a uniform distribution on a unit simplex in dim dimensions. */
    inline std::vector<double> randomSimplexPoint(size_t dim);

    /* Generate n reference points on the unit simplex in dim dimensions (for the NSGA-III algorithm). */
    inline std::vector<std::vector<double>> generateRefPoints(size_t n, size_t dim);

} // namespace genetic_algorithm::detail


/* IMPLEMENTATION */

#include <algorithm>
#include <execution>
#include <random>
#include <limits>
#include <cmath>
#include <cassert>

#include "mo_detail.h"

namespace genetic_algorithm::detail
{
    std::vector<double> randomSimplexPoint(size_t dim)
    {
        assert(dim > 0);

        static thread_local std::minstd_rand0 engine{ std::random_device{}() };
        std::uniform_real_distribution<double> distribution{ 0.0, 1.0 };

        std::vector<double> point;
        point.reserve(dim);

        double sum = 0.0;
        for (size_t i = 0; i < dim; i++)
        {
            point.push_back(-std::log(distribution(engine)));
            sum += point.back();
        }
        for (size_t i = 0; i < dim; i++)
        {
            point[i] /= sum;
        }

        return point;
    }

    std::vector<std::vector<double>> generateRefPoints(size_t n, size_t dim)
    {
        using namespace std;
        assert(n > 0);
        assert(dim > 1);

        /* Generate reference point candidates randomly. */
        size_t k = max(size_t{ 10 }, 2 * dim);
        vector<vector<double>> candidates(k * n - 1);
        generate(candidates.begin(), candidates.end(), [&dim]() { return randomSimplexPoint(dim); });

        vector<vector<double>> refs;
        refs.reserve(n);

        /* The first ref point can be random. */
        refs.push_back(randomSimplexPoint(dim));

        vector<double> min_distances(candidates.size(), numeric_limits<double>::infinity());
        while (refs.size() < n)
        {
            /* Calc the distance of each candidate to the closest ref point. */
            transform(execution::par_unseq, candidates.begin(), candidates.end(), min_distances.begin(), min_distances.begin(),
            [&refs](const vector<double>& candidate, double dmin)
            {
                double d = euclideanDistanceSq(candidate, refs.back());
                return min(dmin, d);
            });

            /* Add the candidate with highest min_distance to the refs. */
            size_t argmax = static_cast<size_t>(max_element(min_distances.begin(), min_distances.end()) - min_distances.begin());
            refs.push_back(move(candidates[argmax]));

            /* Remove the added candidate and the corresponding min_distance. */
            swap(candidates[argmax], candidates.back());
            candidates.pop_back();
            swap(min_distances[argmax], min_distances.back());
            min_distances.pop_back();
        }

        return refs;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_REFERENCE_POINTS_H