/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ALGORITHM_ANY_OBJECTIVE_IMPL_HPP
#define GAPP_ALGORITHM_ANY_OBJECTIVE_IMPL_HPP

#include "any_objective.decl.hpp"
#include "../core/ga_info.hpp"
#include "../core/population.hpp"
#include "../utility/utility.hpp"

namespace gapp::algorithm
{
    template<typename SOA, typename MOA>
    void AnyObjective<SOA, MOA>::initializeImpl(const GaInfo& ga)
    {
        GAPP_ASSERT(ga.num_objectives() > 0);

        if (ga.num_objectives() == 1)
        {
            algorithm_ = SOA{};
            static_cast<Algorithm*>(std::get_if<SOA>(&algorithm_))->initializeImpl(ga);
        }
        else
        {
            algorithm_ = MOA{};
            static_cast<Algorithm*>(std::get_if<MOA>(&algorithm_))->initializeImpl(ga);
        }
    }

    template<typename SOA, typename MOA>
    void AnyObjective<SOA, MOA>::prepareSelectionsImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(!algorithm_.valueless_by_exception());

        if (Algorithm* algorithm = std::get_if<SOA>(&algorithm_))
        {
            algorithm->prepareSelectionsImpl(ga, pop);
            return;
        }
        if (Algorithm* algorithm = std::get_if<MOA>(&algorithm_))
        {
            algorithm->prepareSelectionsImpl(ga, pop);
            return;
        }

        GAPP_UNREACHABLE();
    }

    template<typename SOA, typename MOA>
    size_t AnyObjective<SOA, MOA>::selectImpl(const GaInfo& ga, const PopulationView& pop) const
    {
        GAPP_ASSERT(!algorithm_.valueless_by_exception());

        if (const Algorithm* algorithm = std::get_if<SOA>(&algorithm_))
        {
            return algorithm->selectImpl(ga, pop);
        }
        if (const Algorithm* algorithm = std::get_if<MOA>(&algorithm_))
        {
            return algorithm->selectImpl(ga, pop);
        }

        GAPP_UNREACHABLE();
    }

    template<typename SOA, typename MOA>
    small_vector<size_t> AnyObjective<SOA, MOA>::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(!algorithm_.valueless_by_exception());

        if (Algorithm* algorithm = std::get_if<SOA>(&algorithm_))
        {
            return algorithm->nextPopulationImpl(ga, pop);
        }
        if (Algorithm* algorithm = std::get_if<MOA>(&algorithm_))
        {
            return algorithm->nextPopulationImpl(ga, pop);
        }

        GAPP_UNREACHABLE();
    }

    template<typename SOA, typename MOA>
    small_vector<size_t> AnyObjective<SOA, MOA>::optimalSolutionsImpl(const GaInfo& ga, const PopulationView& pop) const
    {
        GAPP_ASSERT(!algorithm_.valueless_by_exception());

        if (const Algorithm* algorithm = std::get_if<SOA>(&algorithm_))
        {
            return algorithm->optimalSolutionsImpl(ga, pop);
        }
        if (const Algorithm* algorithm = std::get_if<MOA>(&algorithm_))
        {
            return algorithm->optimalSolutionsImpl(ga, pop);
        }

        GAPP_UNREACHABLE();
    }

} // namespace gapp::algorithm

#endif // !GAPP_ALGORITHM_ANY_OBJECTIVE_IMPL_HPP
