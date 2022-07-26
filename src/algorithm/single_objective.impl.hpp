/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_IMPL_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_IMPL_HPP

#include "single_objective.decl.hpp"
#include "../core/ga_info.hpp"
#include <algorithm>
#include <utility>
#include <cassert>

namespace genetic_algorithm::algorithm
{
    template<selection::SelectionType S, update::UpdaterType U>
    SingleObjective<S, U>::SingleObjective(S selection, U updater) :
        selection_(std::move(selection)),
        updater_(std::move(updater))
    {
    }

    template<selection::SelectionType S, update::UpdaterType U>
    void SingleObjective<S, U>::initialize(const GaInfo& ga)
    {
        assert(ga.num_objectives() == 1);

        selection_.initialize(ga);
    }

    template<selection::SelectionType S, update::UpdaterType U>
    void SingleObjective<S, U>::prepareSelections(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        selection_.prepareSelections(ga, fmat);
    }

    template<selection::SelectionType S, update::UpdaterType U>
    size_t SingleObjective<S, U>::select(const GaInfo& ga, const FitnessMatrix& fmat) const
    {
        return selection_.select(ga, fmat);
    }

    template<selection::SelectionType S, update::UpdaterType U>
    std::vector<size_t> SingleObjective<S, U>::nextPopulation(const GaInfo& ga,
                                                              FitnessMatrix::const_iterator first,
                                                              FitnessMatrix::const_iterator children_first,
                                                              FitnessMatrix::const_iterator last)
    {
        assert(ga.num_objectives() == 1);
        assert(std::all_of(first, last, [](const FitnessVector& f) { return f.size() == 1; }));

        return updater_(ga, first, children_first, last);
    }

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_IMPL_HPP