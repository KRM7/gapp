/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_MUTATION_MIXED_DECL_HPP
#define GAPP_MUTATION_MIXED_DECL_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/type_traits.hpp"
#include <array>
#include <tuple>
#include <type_traits>
#include <cstddef>

namespace gapp::mutation
{
    /**
    * The mixed mutation template that is used as the mutation operator in the
    * MixedGene GAs.
    *
    * The mixed mutation consists of a separate component mutation for each of the
    * gene types in the mixed gene. These component mutations are applied separately
    * to the appropriate chromosomes of the mixed gene candidates when performing the
    * mutations.
    *
    * The component mutations are independent of each other, and each of them must be
    * a valid mutation operator that could be used for the given gene type.
    *
    * @tparam Ts The types of the component mutations. Must consist of mutations for
    *   unique gene types. None of the mutation types should be a mutation for some
    *   MixedGene<> type.
    */
    template<typename... Ts>
    class Mixed final : public Mutation<MixedGene<typename Ts::GeneType...>>
    {
    public:
        static_assert((detail::is_derived_from_spec_of_v<Ts, Mutation> && ...),
          "The components of the mixed mutation operator must be valid mutations.");

        using Base = Mutation<MixedGene<typename Ts::GeneType...>>;

        /** The number of component mutations the mixed mutation is composed of. */
        static constexpr size_t N = Base::N;

        /**
        * Create a mixed mutation operator by default constructing all of the component
        * mutations.
        */
        constexpr explicit Mixed() = default;

        /**
        * Create a mixed mutation operator from the specified component mutations.
        * The order of the component mutations must match the order of the gene types in
        * the mixed gene type that the mixed mutation operator is going to be used for.
        *
        * @param mutations The component mutations that will comprise the mixed mutation.
        */
        constexpr explicit Mixed(Ts... mutations);

        /**
        * Set the mutation probability used for each of the component mutations to the
        * same value.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void mutation_rates(Probability pm) noexcept override;

        /**
        * Set the mutation probability used for each of the component mutations
        * individually. The order of the probabilities should match the order of the
        * component mutations.
        *
        * @param pms The mutation probabilities. They must all be in the closed interval [0.0, 1.0].
        */
        constexpr void mutation_rates(const std::array<Probability, N>& pms) noexcept override;

        /**
        * @returns The mutation rates set for the component mutations. The order of the
        *   probabilities in the returned array will match the order of the component mutations.
        */
        [[nodiscard]]
        constexpr std::array<Probability, N> mutation_rates() const noexcept override;

    private:
        using typename Base::GeneType;

        void initialize(const GaInfo& ga) override;
        void mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const override;

        constexpr void mutation_rate_impl(Probability pc, size_t idx) noexcept override;
        constexpr Probability mutation_rate_impl(size_t idx) const noexcept override;

        constexpr bool allow_variable_length_impl(size_t idx) const noexcept override;

        void* component_impl(size_t idx) noexcept override;
        const void* component_impl(size_t idx) const noexcept override;

        std::tuple<Ts...> components_;
    };

} // namespace gapp::mutation

#endif // !GAPP_MUTATION_MIXED_DECL_HPP
