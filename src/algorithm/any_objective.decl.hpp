/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ALGORITHM_ANY_OBJECTIVE_DECL_HPP
#define GAPP_ALGORITHM_ANY_OBJECTIVE_DECL_HPP

#include "algorithm_base.hpp"
#include "single_objective.hpp"
#include "nsga3.hpp"
#include <variant>
#include <type_traits>

namespace gapp::algorithm
{
    /**
    * A simple wrapper algorithm for a single- and a multi-objective algorithm.
    * The single objective algorithm will be used for single-objective problems,
    * while the multi-objective one will be used for multi-objective ones.
    * 
    * This is intended to help turning algorithms which only work for single- or
    * multi-objective problems into algorithms which can be used for any problem
    * type.
    * 
    * @tparam SOA The algorithm type to use for single-objective problems.
    *   Must be default constructible.
    * @tparam MOA The algorithm type to use for multi-objective problems.
    *   Must be default constructible.
    */
    template<typename SOA = algorithm::SingleObjective, typename MOA = algorithm::NSGA3>
    class AnyObjective final : public Algorithm
    {
        static_assert(std::is_default_constructible_v<SOA> && std::is_default_constructible_v<MOA>);

        void initializeImpl(const GaInfo& ga) override;

        void prepareSelectionsImpl(const GaInfo& ga, const PopulationView& pop) override;
        size_t selectImpl(const GaInfo& ga, const PopulationView& pop) const override;

        small_vector<size_t> nextPopulationImpl(const GaInfo& ga, const PopulationView& pop) override;
        small_vector<size_t> optimalSolutionsImpl(const GaInfo& ga, const PopulationView& pop) const override;

        std::variant<SOA, MOA> algorithm_;
    };

} // namespace gapp::algorithm

#endif // !GAPP_ALGORITHM_ANY_OBJECTIVE_DECL_HPP
