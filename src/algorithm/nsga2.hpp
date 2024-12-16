/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ALGORITHM_NSGA2_HPP
#define GAPP_ALGORITHM_NSGA2_HPP

#include "algorithm_base.hpp"
#include "../utility/small_vector.hpp"
#include <vector>
#include <cstddef>

namespace gapp::algorithm
{
    /**
    * Non-dominated sorting genetic algorithm (NSGA-II), used for multi-objective optimization.
    * This algorithm doesn't work for single-objective problems.
    * 
    * The aim of the algorithm is to find a set of solutions which are well spread out
    * along the entire pareto-front in the objective-space.
    * 
    * The algorithm uses a non-dominated sorting method to sort the candidates of the population
    * into distinct pareto fronts, and then selects the candidates of the best fronts for the
    * population of the next generation.
    * Candidates that belong to the same front are ranked according to their crowding distances,
    * which are their distances from the neigbouring solutions in the objective space.
    * 
    * The algorithm uses a selection operator that selects solutions for the crossovers
    * based on these same criteria (their pareto ranks and crowding distances).
    * 
    * The algorithm assumes fitness maximization, and has no parameters.
    * 
    * @see
    *  Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II."
    *  IEEE transactions on evolutionary computation 6, no. 2 (2002): 182-197.
    */
    class NSGA2 final : public Algorithm
    {
    private:
        void initializeImpl(const GaInfo& ga) override;
        void prepareSelectionsImpl(const GaInfo&, const PopulationView&) override {}
        size_t selectImpl(const GaInfo& ga, const PopulationView& pop) const override;

        small_vector<size_t> nextPopulationImpl(const GaInfo& ga, const PopulationView& pop) override;

        small_vector<size_t> optimalSolutionsImpl(const GaInfo& ga, const PopulationView& pop) const override;

        std::vector<size_t> ranks_;
        std::vector<double> dists_;
    };

} // namespace gapp::algorithm

#endif // !GAPP_ALGORITHM_NSGA2_HPP
