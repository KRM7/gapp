/* Copyright (c) 2022 Kriszti�n Rug�si. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_HPP

#include "algorithm_base.hpp"
#include "soga_selection.hpp"
#include "pop_update.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include <utility>
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    /**
    * A generic algorithm for single-objective optimization. \n
    * The algorithm combines a selection method and a population update
    * method. The selection method is used to select candidates from the populations
    * for crossover, while the population update method is used to create the population
    * for the next generation of the algorithm from the combined parent and child populations
    * of the current generation.
    * 
    * @tparam Selection The type of the selection method used.
    * @tparam Updater The type of the population update method used.
    */
    template<selection_::Selection Selection, pop_update::Updater PopUpdater>
    class SingleObjective final : public Algorithm
    {
    public:
        /**
        * Create a single objective algorithm.
        * 
        * @param selection The selection method to use in the algorithm.
        * @param updater The method used to update the populations between generations of the algorithm.
        */
        SingleObjective(Selection selection = selection_::Tournament{},
                        PopUpdater updater  = pop_update::KeepBest{})
            : selection_(std::move(selection)), updater_(std::move(updater))
        {}

        /** @returns The selection method of the algorithm. */
        [[nodiscard]]
        Selection& selection_method() const noexcept { return selection_; }

        /** @returns The population update method of the algorithm. */
        [[nodiscard]]
        PopUpdater& update_method() const noexcept { return updater_; }


        void initialize(const GaInfo& ga) override
        {
            selection_.initialize(ga);
        }

        void prepareSelections(const GaInfo& ga, const FitnessMatrix& population_fmat) override
        {
            selection_.prepareSelections(ga, population_fmat);
        }

        size_t select(const GaInfo& ga, const FitnessMatrix& population_fmat) override
        {
            return selection_.select(ga, population_fmat);
        }

        std::vector<size_t> nextPopulation(const GaInfo& ga,
                                           FitnessMatrix::const_iterator first,
                                           FitnessMatrix::const_iterator children_first,
                                           FitnessMatrix::const_iterator last) override
        {
            return updater_.nextPopulation(ga, population_fmat);
        }

    private:
        Selection selection_;
        PopUpdater updater_;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_HPP