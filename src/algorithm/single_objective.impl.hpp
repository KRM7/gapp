/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_IMPL_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_IMPL_HPP

#include "single_objective.decl.hpp"
#include "../core/ga_info.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <utility>
#include <stdexcept>

namespace genetic_algorithm::algorithm
{
    template<selection::SelectionType S, update::UpdaterType U>
    SingleObjective<S, U>::SingleObjective(S selection, U updater) noexcept :
        selection_(std::move(selection)),
        updater_(std::move(updater))
    {
    }

    template<selection::SelectionType S, update::UpdaterType U>
    void SingleObjective<S, U>::initializeImpl(const GaInfo& ga)
    {
        if (ga.num_objectives() != 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be 1 for the single-objective algorithms.");
        }

        selection_.initialize(ga);
    }

    template<selection::SelectionType S, update::UpdaterType U>
    void SingleObjective<S, U>::prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        selection_.prepareSelections(ga, fmat);
    }

    template<selection::SelectionType S, update::UpdaterType U>
    size_t SingleObjective<S, U>::selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const
    {
        return selection_.select(ga, fmat);
    }

    template<selection::SelectionType S, update::UpdaterType U>
    std::vector<size_t> SingleObjective<S, U>::nextPopulationImpl(const GaInfo& ga,
                                                                  FitnessMatrix::const_iterator first,
                                                                  FitnessMatrix::const_iterator children_first,
                                                                  FitnessMatrix::const_iterator last)
    {
        assert(size_t(children_first - first) == ga.population_size());
        assert(size_t(last - children_first) >= ga.population_size());
        assert(std::all_of(first, last, [](const FitnessVector& fvec) { return fvec.size() == 1; }));

        if (ga.num_objectives() != 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be 1 for the single-objective algorithms.");
        }

        return updater_(ga, first, children_first, last);
    }

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_IMPL_HPP