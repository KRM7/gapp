/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_BASE_IMPL_HPP
#define GAPP_CROSSOVER_BASE_IMPL_HPP

#include "crossover_base.decl.hpp"
#include "../core/ga_info.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/rng.hpp"
#include "../utility/type_traits.hpp"
#include "../utility/utility.hpp"
#include <type_traits>
#include <utility>

namespace gapp::crossover
{
    template<typename T>
    CandidatePair<T> Crossover<T>::operator()(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        GAPP_ASSERT(parent1.is_evaluated() && parent2.is_evaluated());
        GAPP_ASSERT(parent1.fitness.size() == ga.num_objectives());
        GAPP_ASSERT(parent2.fitness.size() == ga.num_objectives());
        GAPP_ASSERT(allow_variable_chrom_length() || parent1.chrom_len() == parent2.chrom_len());

        if constexpr (is_bounded_gene_v<T>)
        {
            GAPP_ASSERT(parent1.chromosome.size() == parent1.gene_bounds.size());
            GAPP_ASSERT(detail::equal(parent1.gene_bounds, parent2.gene_bounds));
        }

        /* Only need to perform the crossover with the set pc probability. Return early with (1 - pc) probability. */
        if (rng::randomReal() >= pc_)
        {
            return { parent1, parent2 };
        }

        /*
        * If the parents are the same, the crossover doesn't need to be performed.
        * This assumes that with 2 parents that have the same chromosomes, the children would be the same
        * as the parents, which is true for every crossover operator implemented, but
        * could be an issue for user defined crossovers.
        */
        if (parent1 == parent2)
        {
            return { parent1, parent2 };
        }

        /* Perform the actual crossover. */
        auto [child1, child2] = crossover(ga, parent1, parent2);

        GAPP_ASSERT(allow_variable_chrom_length() || child1.chromosome.size() == parent1.chromosome.size(),
                  "The crossover created a candidate with incorrect chromosome length.");
        GAPP_ASSERT(allow_variable_chrom_length() || child2.chromosome.size() == parent2.chromosome.size(),
                  "The crossover created a candidate with incorrect chromosome length.");

        if constexpr (is_bounded_gene_v<T>)
        {
            child1.gene_bounds = parent1.gene_bounds;
            child2.gene_bounds = parent2.gene_bounds;
        }

        return { std::move(child1), std::move(child2) };
    }

    template<typename... Ts>
    auto Crossover<MixedGene<Ts...>>::operator()(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const
        -> CandidatePair<GeneType>
    {
        return crossover(ga, parent1, parent2);
    }

    template<typename... Ts>
    template<typename GeneType>
    constexpr void Crossover<MixedGene<Ts...>>::crossover_rate(Probability pc) noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        crossover_rate_impl(pc, detail::index_of_type_v<GeneType, Ts...>);
    }

    template<typename... Ts>
    template<typename GeneType>
    constexpr Probability Crossover<MixedGene<Ts...>>::crossover_rate() const noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...), 
          "The GeneType must be one of the gene types in the mixed gene.");

        return crossover_rate_impl(detail::index_of_type_v<GeneType, Ts...>);
    }

    template<typename... Ts>
    template<typename GeneType>
    constexpr bool Crossover<MixedGene<Ts...>>::allow_variable_chrom_length() const noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return allow_variable_length_impl(detail::index_of_type_v<GeneType, Ts...>);
    }

    template<typename... Ts>
    template<typename GeneType>
    Crossover<GeneType>& Crossover<MixedGene<Ts...>>::component() & noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return *static_cast<Crossover<GeneType>*>(component_impl(detail::index_of_type_v<GeneType, Ts...>));
    }

    template<typename... Ts>
    template<typename GeneType>
    const Crossover<GeneType>& Crossover<MixedGene<Ts...>>::component() const& noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return *static_cast<const Crossover<GeneType>*>(component_impl(detail::index_of_type_v<GeneType, Ts...>));
    }

} // namespace gapp::crossover

#endif // !GAPP_CROSSOVER_BASE_IMPL_HPP
