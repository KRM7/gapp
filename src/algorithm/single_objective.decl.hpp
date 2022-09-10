/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_DECL_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_DECL_HPP

#include "algorithm_base.hpp"
#include "single_objective.fwd.hpp"
#include "soga_selection.hpp"
#include "pop_update.hpp"

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
    template<selection::SelectionType Selection = selection::Tournament, update::UpdaterType Updater = update::KeepBest>
    class SingleObjective final : public Algorithm
    {
    public:
        /**
        * Create a single objective algorithm.
        * 
        * @param selection The selection method to use in the algorithm.
        * @param updater The method used to update the populations between generations of the algorithm.
        */
        SingleObjective(Selection selection = Selection{}, Updater updater = Updater{}) noexcept;

        /** @returns The selection method of the algorithm. */
        [[nodiscard]]
        Selection& selection_method() const noexcept { return selection_; }

        /** @returns The population update method of the algorithm. */
        [[nodiscard]]
        Updater& update_method() const noexcept { return updater_; }

    private:
        Selection selection_;
        Updater updater_;

        void initializeImpl(const GaInfo& ga) override;
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<size_t> nextPopulationImpl(const GaInfo& ga,
                                               FitnessMatrix::const_iterator first,
                                               FitnessMatrix::const_iterator children_first,
                                               FitnessMatrix::const_iterator last) override;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_DECL_HPP