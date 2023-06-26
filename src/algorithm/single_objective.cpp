/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "single_objective.hpp"
#include "soga_selection.hpp"
#include "soga_replacement.hpp"
#include "../core/ga_info.hpp"
#include "../utility/functional.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <memory>
#include <utility>

namespace gapp::algorithm
{
    SingleObjective::SingleObjective(std::unique_ptr<selection::Selection> selection) :
        SingleObjective(std::move(selection), std::make_unique<DefaultReplacement>())
    {}

    SingleObjective::SingleObjective(std::unique_ptr<selection::Selection> selection, std::unique_ptr<replacement::Replacement> replacement) :
        selection_(std::move(selection)), replacement_(std::move(replacement))
    {
        GAPP_ASSERT(selection_, "The selection method can't be a nullptr.");
        GAPP_ASSERT(replacement_, "The population replacement method can't be a nullptr.");
    }

    SingleObjective::SingleObjective(SelectionCallable selection) :
        SingleObjective(std::make_unique<selection::Lambda>(std::move(selection)))
    {}

    SingleObjective::SingleObjective(SelectionCallable selection, ReplacementCallable replacement) :
        SingleObjective(std::make_unique<selection::Lambda>(std::move(selection)), std::make_unique<replacement::Lambda>(std::move(replacement)))
    {}

    void SingleObjective::selection_method(std::unique_ptr<selection::Selection> selection)
    {
        GAPP_ASSERT(selection_, "The selection method can't be a nullptr.");

        selection_ = std::move(selection);
    }

    void SingleObjective::selection_method(SelectionCallable f)
    {
        selection_ = std::make_unique<selection::Lambda>(std::move(f));
    }

    void SingleObjective::replacement_method(std::unique_ptr<replacement::Replacement> replacement)
    {
        GAPP_ASSERT(replacement_, "The population replacement method can't be a nullptr.");

        replacement_ = std::move(replacement);
    }

    void SingleObjective::replacement_method(ReplacementCallable f)
    {
        replacement_ = std::make_unique<replacement::Lambda>(std::move(f));
    }

    void SingleObjective::initializeImpl(const GaInfo& ga)
    {
        GAPP_ASSERT(selection_);
        GAPP_ASSERT(ga.num_objectives() == 1, "The number of objectives must be 1 for the single-objective algorithms.");

        selection_->initializeImpl(ga);
    }

    std::vector<size_t> SingleObjective::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last)
    {
        GAPP_ASSERT(replacement_);
        GAPP_ASSERT(ga.num_objectives() == 1, "The number of objectives must be 1 for the single-objective algorithms.");

        return replacement_->nextPopulationImpl(ga, first, children_first, last);
    }

} // namespace gapp::algorithm