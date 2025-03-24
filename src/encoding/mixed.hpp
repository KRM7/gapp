/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ENCODING_MIXED_HPP
#define GAPP_ENCODING_MIXED_HPP

#include "gene_types.hpp"
#include "../core/ga_traits.hpp"
#include "../core/ga_base.hpp"
#include "../core/candidate.hpp"
#include "../crossover/mixed.hpp"
#include "../mutation/mixed.hpp"

namespace gapp
{
    template<typename... Ts>
    struct GaTraits<MixedGene<Ts...>>
    {
        using DefaultCrossover = crossover::Mixed<typename GaTraits<Ts>::DefaultCrossover...>;
        using DefaultMutation  = mutation::Mixed<typename GaTraits<Ts>::DefaultMutation...>;
    };

    template<typename... Ts>
    class MixedGA final : public GA<MixedGene<Ts...>>
    {
    public:
        using GA<MixedGene<Ts...>>::GA;
    };

} // namespace gapp

#endif // !GAPP_ENCODING_MIXED_HPP
