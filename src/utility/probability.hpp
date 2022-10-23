/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_PROBABILITY_HPP
#define GA_UTILITY_PROBABILITY_HPP

#include "utility.hpp"
#include <stdexcept>

namespace genetic_algorithm
{
    /* Class representing a probability value in the closed range [0.0, 1.0]. */
    class Probability
    {
    public:
        constexpr /* implicit */ Probability(double p)
        {
            if (!(0.0 <= p && p <= 1.0))
            {
                GA_THROW(std::invalid_argument, "Probabilities must be in the closed range [0.0, 1.0].");
            }
            p_ = p;
        }

        constexpr /* implicit */ operator double() const noexcept { return p_; }

    private:
        double p_;
    };

} // namespace genetic_algorithm

#endif // !GA_UTILITY_PROBABILITY_HPP