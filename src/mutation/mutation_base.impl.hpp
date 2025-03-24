/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_MUTATION_BASE_IMPL_HPP
#define GAPP_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../core/ga_info.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/type_traits.hpp"
#include "../utility/utility.hpp"
#include <type_traits>
#include <utility>

namespace gapp::mutation
{
    template<typename T>
    void Mutation<T>::operator()(const GaInfo& ga, Candidate<T>& candidate) const
    {
        GAPP_ASSERT(!candidate.is_evaluated() || candidate.fitness.size() == ga.num_objectives());
        GAPP_ASSERT(allow_variable_chrom_length() || candidate.chromosome.size() == ga.chrom_len<T>());

        mutate(ga, candidate, candidate.chromosome);

        GAPP_ASSERT(allow_variable_chrom_length() || candidate.chromosome.size() == ga.chrom_len<T>(),
          "The mutation resulted in a candidate with incorrect chromosome length.");
    }

    template<typename... Ts>
    void Mutation<MixedGene<Ts...>>::operator()(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        mutate(ga, candidate);
    }

    template<typename... Ts>
    template<typename GeneType>
    constexpr void Mutation<MixedGene<Ts...>>::mutation_rate(Probability pm) noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        mutation_rate_impl(pm, detail::index_of_type_v<GeneType, Ts...>);
    }

    template<typename... Ts>
    template<typename GeneType>
    constexpr Probability Mutation<MixedGene<Ts...>>::mutation_rate() const noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return mutation_rate_impl(detail::index_of_type_v<GeneType, Ts...>);
    }

    template<typename... Ts>
    template<typename GeneType>
    constexpr bool Mutation<MixedGene<Ts...>>::allow_variable_chrom_length() const noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return allow_variable_length_impl(detail::index_of_type_v<GeneType, Ts...>);
    }

    template<typename... Ts>
    template<typename GeneType>
    Mutation<GeneType>& Mutation<MixedGene<Ts...>>::component() & noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return *static_cast<Mutation<GeneType>*>(component_impl(detail::index_of_type_v<GeneType, Ts...>));
    }

    template<typename... Ts>
    template<typename GeneType>
    const Mutation<GeneType>& Mutation<MixedGene<Ts...>>::component() const& noexcept
    {
        static_assert((std::is_same_v<GeneType, Ts> || ...),
          "The GeneType must be one of the gene types in the mixed gene.");

        return *static_cast<const Mutation<GeneType>*>(component_impl(detail::index_of_type_v<GeneType, Ts...>));
    }

} // namespace gapp::mutation

#endif // !GAPP_MUTATION_BASE_IMPL_HPP
