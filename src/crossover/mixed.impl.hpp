/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_MIXED_IMPL_HPP
#define GAPP_CROSSOVER_MIXED_IMPL_HPP

#include "crossover_base.hpp"
#include "../core/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/type_traits.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <array>
#include <tuple>
#include <type_traits>
#include <memory>
#include <utility>

namespace gapp::crossover
{
    template<typename... Ts>
    constexpr Mixed<Ts...>::Mixed(Ts... crossovers) :
        components_(std::move(crossovers)...)
    {}

    template<typename... Ts>
    constexpr void Mixed<Ts...>::crossover_rate_impl(Probability pc, size_t idx) noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        ((detail::index_of_type_v<Ts, Ts...> == idx ? std::get<Ts>(components_).crossover_rate(pc) : void()), ...);
    }

    template<typename... Ts>
    constexpr void Mixed<Ts...>::crossover_rates(Probability pc) noexcept
    {
        (std::get<Ts>(components_).crossover_rate(pc), ...);
    }

    template<typename... Ts>
    constexpr void Mixed<Ts...>::crossover_rates(const std::array<Probability, N>& pcs) noexcept
    {
        (std::get<Ts>(components_).crossover_rate(pcs[detail::index_of_type_v<Ts, Ts...>]), ...);
    }

    template<typename... Ts>
    constexpr Probability Mixed<Ts...>::crossover_rate_impl(size_t idx) const noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        return crossover_rates()[idx];
    }

    template<typename... Ts>
    constexpr auto Mixed<Ts...>::crossover_rates() const noexcept -> std::array<Probability, N>
    {
        return { std::get<Ts>(components_).crossover_rate()... };
    }

    template<typename... Ts>
    constexpr bool Mixed<Ts...>::allow_variable_length_impl(size_t idx) const noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        return std::array{ std::get<Ts>(components_).allow_variable_chrom_length()... }[idx];
    }

    template<typename... Ts>
    void* Mixed<Ts...>::component_impl(size_t idx) noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        return std::array{ static_cast<void*>(static_cast<Crossover<typename Ts::GeneType>*>(std::addressof(std::get<Ts>(components_))))... }[idx];
    }

    template<typename... Ts>
    const void* Mixed<Ts...>::component_impl(size_t idx) const noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        return std::array{ static_cast<const void*>(static_cast<const Crossover<typename Ts::GeneType>*>(std::addressof(std::get<Ts>(components_))))... }[idx];
    }

    template<typename... Ts>
    void Mixed<Ts...>::initialize(const GaInfo& ga)
    {
        (static_cast<Crossover<typename Ts::GeneType>&>(std::get<Ts>(components_)).initialize(ga), ...);
    }

    template<typename... Ts>
    auto Mixed<Ts...>::crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const
        -> CandidatePair<GeneType>
    {
        std::tuple children{ std::get<Ts>(components_)(ga, parent1, parent2)... };

        Candidate<GeneType> child1({ std::move(std::get<CandidatePair<typename Ts::GeneType>>(children).first)... });
        Candidate<GeneType> child2({ std::move(std::get<CandidatePair<typename Ts::GeneType>>(children).second)... });

        return { std::move(child1), std::move(child2) };
    }

} // namespace gapp::crossover

#endif // !GAPP_CROSSOVER_MIXED_IMPL_HPP
