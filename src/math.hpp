#ifndef GA_MATH_HPP
#define GA_MATH_HPP

#include "utils.h"

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