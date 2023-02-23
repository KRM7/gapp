/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "single_objective.hpp"
#include "soga_selection.hpp"
#include "soga_update.hpp"
#include "../core/ga_info.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>

namespace genetic_algorithm::algorithm
{
    SingleObjective::SingleObjective(SelectionCallable selection) :
        SingleObjective(std::make_unique<selection::Lambda>(std::move(selection)))
    {}

    SingleObjective::SingleObjective(SelectionCallable selection, UpdateCallable updater) :
        SingleObjective(std::make_unique<selection::Lambda>(std::move(selection)), std::make_unique<update::Lambda>(std::move(updater)))
    {}

    void SingleObjective::selection_method(SelectionCallable f)
    {
        selection_ = std::make_unique<selection::Lambda>(std::move(f));
    }

    void SingleObjective::update_method(UpdateCallable f)
    {
        updater_ = std::make_unique<update::Lambda>(std::move(f));
    }

    void SingleObjective::initializeImpl(const GaInfo& ga)
    {
        if (ga.num_objectives() != 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be 1 for the single-objective algorithms.");
        }

        selection_->initializeImpl(ga);
    }

    std::vector<size_t> SingleObjective::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last)
    {
        assert(std::all_of(first, last, detail::is_size(1)));

        if (ga.num_objectives() != 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be 1 for the single-objective algorithms.");
        }

        return updater_->nextPopulationImpl(ga, first, children_first, last);
    }

} // namespace genetic_algorithm::algorithm