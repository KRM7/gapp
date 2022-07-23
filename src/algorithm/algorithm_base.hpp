/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_ALGORITHM_BASE_HPP
#define GA_ALGORITHM_ALGORITHM_BASE_HPP

#include "../population/population.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;
}

/** Algorithm types that can be used in the genetic algorithms. */
namespace genetic_algorithm::algorithm
{
    /** ... */
    class Algorithm
    {
    public:
        using FitnessVector = detail::FitnessVector;
        using FitnessMatrix = detail::FitnessMatrix;

        /** ... */
        virtual void initialize(const GaInfo& ga) = 0;

        /** ... */
        virtual void prepareSelections(const GaInfo& ga, const FitnessMatrix& population_fmat) = 0;

        /** ... */
        virtual size_t select(const GaInfo& ga, const FitnessMatrix& population_fmat) = 0;

        /** ... */
        virtual std::vector<size_t> nextPopulation(const GaInfo& ga,
                                                   FitnessMatrix::const_iterator first,
                                                   FitnessMatrix::const_iterator children_first,
                                                   FitnessMatrix::const_iterator last) = 0;

        Algorithm()                             = default;
        Algorithm(const Algorithm&)             = default;
        Algorithm(Algorithm&&)                  = default;
        Algorithm& operator=(const Algorithm&)  = default;
        Algorithm& operator=(Algorithm&&)       = default;
        virtual ~Algorithm()                    = default;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_ALGORITHM_BASE_HPP