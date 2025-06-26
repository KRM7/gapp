/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_MIXED_DECL_HPP
#define GAPP_CROSSOVER_MIXED_DECL_HPP

#include "crossover_base.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/type_traits.hpp"
#include "../utility/bounded_value.hpp"
#include <array>
#include <tuple>
#include <type_traits>
#include <cstddef>

namespace gapp::crossover
{
    /**
    * The mixed crossover template that is used as the crossover operator in the
    * MixedGene GAs.
    * 
    * The mixed crossover consists of a separate component crossover for each of the
    * gene types in the mixed gene. These component crossovers are applied separately
    * to the appropriate chromosomes of the mixed gene candidates in order to create
    * the child candidates.
    * 
    * The component crossovers are independent of each other, and each of them must be
    * a valid crossover operator that could be used for the given gene type.
    * 
    * @tparam Ts The types of the component crossovers. Must consist of crossovers for
    *   unique gene types. None of the crossover types should be a crossover for some
    *   MixedGene<> type.
    */
    template<typename... Ts>
    class Mixed final : public Crossover<MixedGene<typename Ts::GeneType...>>
    {
    public:
        static_assert((detail::is_derived_from_spec_of_v<Ts, Crossover> && ...), 
          "The components of the mixed crossover operator must be valid crossovers.");

        using Base = Crossover<MixedGene<typename Ts::GeneType...>>;

        /** The number of component crossovers the mixed crossover is composed of. */
        static constexpr size_t N = Base::N;

        /**
        * Create a mixed crossover operator by default constructing all of the component
        * crossovers.
        */
        constexpr explicit Mixed() = default;

        /**
        * Create a mixed crossover operator from the specified component crossovers.
        * The order of the component crossovers must match the order of the gene types in
        * the mixed gene type that the mixed crossover operator is going to be used for.
        * 
        * @param crossovers The component crossovers that will comprise the mixed crossover.
        */
        constexpr explicit Mixed(Ts... crossovers);

        /**
        * Set the crossover probability used for each of the component crossovers to the
        * same value.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void crossover_rates(Probability pc) noexcept override;

        /**
        * Set the crossover probability used for each of the component crossovers
        * individually. The order of the probabilities should match the order of the
        * component crossovers.
        *
        * @param pcs The crossover probabilities. They must all be in the closed interval [0.0, 1.0].
        */
        constexpr void crossover_rates(const std::array<Probability, N>& pcs) noexcept override;

        /**
        * @returns The crossover rates set for the component crossovers. The order of the
        *   probabilities in the returned array will match the order of the component crossovers.
        */
        [[nodiscard]]
        constexpr std::array<Probability, N> crossover_rates() const noexcept;

    private:
        using typename Base::GeneType;

        void initialize(const GaInfo& ga) override;
        CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override;

        constexpr void crossover_rate_impl(Probability pc, size_t idx) noexcept override;
        constexpr Probability crossover_rate_impl(size_t idx) const noexcept override;

        constexpr bool allow_variable_length_impl(size_t idx) const noexcept override;

        void* component_impl(size_t idx) noexcept override;
        const void* component_impl(size_t idx) const noexcept override;

        std::tuple<Ts...> components_;
    };

} // namespace gapp::crossover

#endif // !GAPP_CROSSOVER_MIXED_DECL_HPP
