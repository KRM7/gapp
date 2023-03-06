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
#include <utility>

namespace genetic_algorithm::algorithm
{
    SingleObjective::SingleObjective(std::unique_ptr<selection::Selection> selection) :
        SingleObjective(std::move(selection), std::make_unique<DefaultUpdater>())
    {}

    SingleObjective::SingleObjective(std::unique_ptr<selection::Selection> selection, std::unique_ptr<update::Updater> updater) :
        selection_(std::move(selection)), updater_(std::move(updater))
    {
        if (!selection_) GA_THROW(std::invalid_argument, "The selection method can't be a nullptr.");
        if (!updater_)   GA_THROW(std::invalid_argument, "The population update method can't be a nullptr.");
    }

    SingleObjective::SingleObjective(SelectionCallable selection) :
        SingleObjective(std::make_unique<selection::Lambda>(std::move(selection)))
    {}

    SingleObjective::SingleObjective(SelectionCallable selection, UpdateCallable updater) :
        SingleObjective(std::make_unique<selection::Lambda>(std::move(selection)), std::make_unique<update::Lambda>(std::move(updater)))
    {}

    void SingleObjective::selection_method(std::unique_ptr<selection::Selection> selection)
    {
        if (!selection_) GA_THROW(std::invalid_argument, "The selection method can't be a nullptr.");

        selection_ = std::move(selection);
    }

    void SingleObjective::selection_method(SelectionCallable f)
    {
        selection_ = std::make_unique<selection::Lambda>(std::move(f));
    }

    void SingleObjective::update_method(std::unique_ptr<update::Updater> updater)
    {
        if (!updater_) GA_THROW(std::invalid_argument, "The population update method can't be a nullptr.");

        updater_ = std::move(updater);
    }

    void SingleObjective::update_method(UpdateCallable f)
    {
        updater_ = std::make_unique<update::Lambda>(std::move(f));
    }

    void SingleObjective::initializeImpl(const GaInfo& ga)
    {
        GA_ASSERT(selection_);

        if (ga.num_objectives() != 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be 1 for the single-objective algorithms.");
        }

        selection_->initializeImpl(ga);
    }

    std::vector<size_t> SingleObjective::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last)
    {
        GA_ASSERT(updater_);

        if (ga.num_objectives() != 1)
        {
            GA_THROW(std::logic_error, "The number of objectives must be 1 for the single-objective algorithms.");
        }

        return updater_->nextPopulationImpl(ga, first, children_first, last);
    }

} // namespace genetic_algorithm::algorithm