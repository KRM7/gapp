/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_MUTATION_MIXED_IMPL_HPP
#define GAPP_MUTATION_MIXED_IMPL_HPP

#include "mutation_base.hpp"
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

namespace gapp::mutation
{
    template<typename... Ts>
    constexpr Mixed<Ts...>::Mixed(Ts... mutations) :
        components_(std::move(mutations)...)
    {}

    template<typename... Ts>
    constexpr void Mixed<Ts...>::mutation_rate_impl(Probability pm, size_t idx) noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        ((detail::index_of_type_v<Ts, Ts...> == idx ? std::get<Ts>(components_).mutation_rate(pm) : void()), ...);
    }

    template<typename... Ts>
    constexpr void Mixed<Ts...>::mutation_rates(Probability pm) noexcept
    {
        (std::get<Ts>(components_).mutation_rate(pm), ...);
    }

    template<typename... Ts>
    constexpr void Mixed<Ts...>::mutation_rates(const std::array<Probability, N>& pms) noexcept
    {
        (std::get<Ts>(components_).mutation_rate(pms[detail::index_of_type_v<Ts, Ts...>]), ...);
    }

    template<typename... Ts>
    constexpr Probability Mixed<Ts...>::mutation_rate_impl(size_t idx) const noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        return mutation_rates()[idx];
    }

    template<typename... Ts>
    constexpr auto Mixed<Ts...>::mutation_rates() const noexcept -> std::array<Probability, N>
    {
        return { std::get<Ts>(components_).mutation_rate()... };
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

        return std::array{ static_cast<void*>(static_cast<Mutation<typename Ts::GeneType>*>(std::addressof(std::get<Ts>(components_))))... }[idx];
    }

    template<typename... Ts>
    const void* Mixed<Ts...>::component_impl(size_t idx) const noexcept
    {
        GAPP_ASSERT(idx < sizeof...(Ts));

        return std::array{ static_cast<const void*>(static_cast<const Mutation<typename Ts::GeneType>*>(std::addressof(std::get<Ts>(components_))))... }[idx];
    }

    template<typename... Ts>
    void Mixed<Ts...>::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        (std::get<Ts>(components_)(ga, candidate), ...);
    }

} // namespace gapp::mutation

#endif // !GAPP_MUTATION_MIXED_IMPL_HPP
