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

/** Math utility functions. */

#ifndef GA_MATH_HPP
#define GA_MATH_HPP

#include "utils.hpp"

#include <vector>
#include <limits>

namespace genetic_algorithm::detail
{
    /* Equality comparison for floating point numbers. Returns true if lhs is approximately equal to rhs. */
    bool floatIsEqual(double lhs, double rhs, double eps = GA_DEFAULT_EPSILON);

    /* Less than comparison for floating point numbers. Returns true if lhs is definitely less than rhs. */
    bool floatIsLess(double lhs, double rhs, double eps = GA_DEFAULT_EPSILON);

    /* Equality comparison for fp vectors. Returns true if the elements of the vectors are approximately equal. */
    bool floatIsEqual(const std::vector<double>& lhs, const std::vector<double>& rhs, double eps = GA_DEFAULT_EPSILON);

    /* Pareto comparison for fp vectors. Returns true if lhs is dominated by rhs (lhs < rhs) assuming maximization. */
    bool paretoCompareLess(const std::vector<double>& lhs, const std::vector<double>& rhs, double eps = GA_DEFAULT_EPSILON);

    /* Calculate the square of the Euclidean distance between the vectors v1 and v2. */
    double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2);

    /* Calculate the square of the perpendicular distance between a line and a point. */
    double perpendicularDistanceSq(const std::vector<double>& line, const std::vector<double>& point);
}

#endif // !GA_MATH_HPP