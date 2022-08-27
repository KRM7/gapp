/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_PROBABILITY
#define GA_UTILITY_PROBABILITY

#include "utility.hpp"
#include <stdexcept>
#include <compare>

namespace genetic_algorithm
{
    /** Class representing a probability value in the closed range [0.0, 1.0]. */
    class Probability
    {
    public:
        /** Create a probability with a value of 0. */
        constexpr Probability() = default;

        /** Create a probability with a value of p. This value must be in the closed range [0.0, 1.0]. */
        constexpr Probability(double p)
        {
            if (!(0.0 <= p && p <= 1.0))
            {
                GA_THROW(std::invalid_argument, "Probabilities must be in the closed range [0.0, 1.0].");
            }
            p_ = p;
        }

        constexpr operator double() const noexcept { return p_; }           /**< Implicit conversion to double. */
        constexpr const double operator*() const noexcept { return p_; }    /**< Convert to double. */

        constexpr friend auto operator<=>(const Probability&, const Probability&) = default;

    private:
        double p_ = 0.0;
    };

} // namespace genetic_algorithm

#endif // !GA_UTILITY_PROBABILITY