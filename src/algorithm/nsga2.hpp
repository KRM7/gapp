/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_NSGA2_HPP
#define GA_ALGORITHM_NSGA2_HPP

#include "algorithm_base.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    /**
    * Non-dominated sorting genetic algorithm (NSGA-II), used for multi-objective optimization.
    * (This algorithm doesn't work for single-objective problems.) \n
    * The aim of the algorithm is to find a set of solutions which are well-spread out
    * along the entire pareto-front (in the objective-space). \n
    * 
    * The algorithm uses non-dominated sorting to sort the candidates into
    * pareto fronts, and then selects the candidates of the best fronts for the population
    * of the next generation. \n
    * Candidates that belong to the same front are ranked according to their crowding distances,
    * which are their distances from the neigbouring solutions in the objective space. \n
    * 
    * The algorithm uses a selection operator that selects solutions based on these same criteria
    * (their pareto ranks and crowding distances). \n
    * 
    * The selection and population update methods of this algorithm can't be changed. \n
    * Has no parameters. \n
    * 
    * See:
    *  Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II."
    *  IEEE transactions on evolutionary computation 6, no. 2 (2002): 182-197.
    */
    class NSGA2 final : public Algorithm
    {
    private:
        void initializeImpl(const GaInfo& ga) override;
        void prepareSelectionsImpl(const GaInfo&, const FitnessMatrix&) override {}
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& pop) const override;

        std::vector<size_t> nextPopulationImpl(const GaInfo& ga,
                                               FitnessMatrix::const_iterator first,
                                               FitnessMatrix::const_iterator children_first,
                                               FitnessMatrix::const_iterator last) override;

        std::vector<size_t> optimalSolutionsImpl(const GaInfo& ga) const override;

        std::vector<size_t> ranks_;
        std::vector<double> dists_;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_NSGA2_HPP