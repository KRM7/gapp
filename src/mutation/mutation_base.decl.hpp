/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_MUTATION_BASE_DECL_HPP
#define GAPP_MUTATION_BASE_DECL_HPP

#include "../core/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"
#include <array>
#include <optional>
#include <utility>
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

namespace gapp::mutation
{
    /**
    * The base class used for the mutation operators of the GAs.
    * %Mutation operators take a candidate solution, and modify it in some way with a
    * given probability. This probability can be interpreted either per-candidate or
    * per-gene depending on how the particular operator is defined.
    * 
    * New mutation operators should be derived from this class, and they should
    * implement the following virtual method:
    * 
    *   - mutate : Perform the mutation on a single candidate's chromosome.
    * 
    * @tparam T The gene type the mutation operator is defined for.
    */
    template<typename T>
    class Mutation
    {
    public:
        /** The gene type the mutation operator is defined for. */
        using GeneType = T;

        /** Create a mutation operator that will use the default mutation probability. */
        constexpr explicit Mutation() noexcept = default;

        /**
        * Create a mutation operator with the specified mutation probability.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Mutation(Probability pm) noexcept :
            pm_(pm)
        {}

        /**
        * Set the mutation rate used by the operator.
        * 
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void mutation_rate(Probability pm) noexcept { pm_ = pm; };

        /** @returns The mutation rate set for the operator. */
        [[nodiscard]]
        constexpr Probability mutation_rate() const noexcept { return pm_.value_or(1.0); }

        /** @returns True if the operator will use the default mutation rate of the GA. */
        constexpr bool use_default_mutation_rate() const noexcept { return !pm_; }

        /**
        * This method specifies whether the mutation operator supports variable
        * chromosome lengths or not. If variable chromosome lengths are supported,
        * the Candidates passed to the mutation operator are allowed to have
        * chromosome lengths that are different from the chromosome length specified
        * for the GA that the operator is used in. Otherwise the chromosome length of
        * the given gene type must be the same for every Candidate.
        *
        * This method will return false by default. If a particular mutation method allows
        * variable chromosome lengths, it should override this method to return true.
        *
        * @returns True if the mutation operator support variable chromosome lengths.
        */
        constexpr virtual bool allow_variable_chrom_length() const noexcept { return false; }

        /**
        * Perform mutation on a candidate using the set mutation probability.
        * Implemented by mutate().
        *
        * @param ga The genetic algorithm the mutation operator is being used in.
        * @param candidate The candidate to mutate.
        */
        void operator()(const GaInfo& ga, Candidate<T>& candidate) const;


        /** Destructor. */
        virtual ~Mutation()                     = default;

    protected:

        Mutation(const Mutation&)               = default;
        Mutation(Mutation&&)                    = default;
        Mutation& operator=(const Mutation&)    = default;
        Mutation& operator=(Mutation&&)         = default;

    private:

        /**
        * The implementation of the mutation operator. Performs the mutation
        * on the given chromosome in-place with the set probability.
        * This function must handle the mutation probability properly as part
        * of its implementation. The mutated chromosome should be valid candidate
        * solution for the given problem and %GA.
        * 
        * This method will be called exactly once for each child solution
        * in every population.
        *
        * The function must be thread-safe if parallel execution is enabled for the
        * GAs (which is true by default).
        * 
        * @param ga The genetic algorithm the mutation operator is being used in.
        * @param candidate The candidate solution that will be mutated.
        * @param chromosome The chromosome to mutate. This is the chromosome of @p candidate.
        */
        virtual void mutate(const GaInfo& ga, const Candidate<T>& candidate, Chromosome<T>& chromosome) const = 0;

        std::optional<Probability> pm_;
    };

    /**
    * The base class used for the mixed mutation operators.
    * This is effectively the same as the primary template without the mutation
    * probability.
    *
    * @tparam Ts The component gene types of the mixed gene type the mutation operator is defined for.
    */
    template<typename... Ts>
    class Mutation<MixedGene<Ts...>>
    {
    public:
        /** The mixed gene type the mutation operator is defined for. */
        using GeneType = MixedGene<Ts...>;

        /** The number of component mutations the mixed mutation is composed of. */
        static constexpr size_t N = sizeof...(Ts);

        /**
        * Set the mutation probability of the component mutation associated with the
        * specified GeneType.
        *
        * @tparam GeneType The gene type of the component mutation to set the probability of.
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        template<typename GeneType>
        constexpr void mutation_rate(Probability pm) noexcept;

        /**
        * Set the mutation probability used for each of the component mutations to the
        * same value.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr virtual void mutation_rates(Probability pm) noexcept = 0;

        /**
        * Set the mutation probability used for each of the component mutations
        * individually. The order of the probabilities should match the order of the
        * component mutations.
        *
        * @param pms The mutation probabilities. They must all be in the closed interval [0.0, 1.0].
        */
        constexpr virtual void mutation_rates(const std::array<Probability, N>& pms) noexcept = 0;

        /**
        * @returns The mutation probability of the component mutation associated with the specified GeneType.
        * @tparam GeneType The gene type of the component mutation to get the probability of.
        */
        template<typename GeneType>
        [[nodiscard]] constexpr Probability mutation_rate() const noexcept;

        /**
        * @returns The mutation rates set for the component mutations. The order of the
        *   probabilities in the returned array will match the order of the component mutations.
        */
        [[nodiscard]]
        constexpr virtual std::array<Probability, N> mutation_rates() const noexcept = 0;

        /**
        * This method specifies whether the mutation operator supports variable
        * chromosome lengths for the specified gene type or not. If variable chromosome
        * lengths are supported, the Candidates passed to the mutation operator are
        * allowed to have chromosome lengths that are different from the chromosome length
        * specified for the GA that the operator is used in. Otherwise the chromosome length
        * of the given gene type must be the same for every Candidate.
        *
        * @returns True if the mutation operator support variable chromosome lengths for
        *   the specified gene type.
        * @tparam The gene type of the chromosome to check.
        */
        template<typename GeneType>
        [[nodiscard]] constexpr bool allow_variable_chrom_length() const noexcept;

        /**
        * @returns The component mutation associated with the specified gene type.
        * @tparam GeneType The gene type of the component mutation to get.
        */
        template<typename GeneType>
        [[nodiscard]] Mutation<GeneType>& component() & noexcept;

        /**
        * @returns The component mutation associated with the specified gene type.
        * @tparam GeneType The gene type of the component mutation to get.
        */
        template<typename GeneType>
        [[nodiscard]] const Mutation<GeneType>& component() const& noexcept;

        /**
        * Perform the mutation on a candidate solution.
        * Implemented by mutate().
        *
        * @param ga The genetic algorithm the mutation operator is being used in.
        * @param candidate The candidate to mutate.
        */
        void operator()(const GaInfo& ga, Candidate<GeneType>& candidate) const;


        /** Destructor. */
        virtual ~Mutation()                  = default;

    protected:

        Mutation()                           = default;
        Mutation(const Mutation&)            = default;
        Mutation(Mutation&&)                 = default;
        Mutation& operator=(const Mutation&) = default;
        Mutation& operator=(Mutation&&)      = default;

    private:
        virtual void mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const = 0;

        constexpr virtual void mutation_rate_impl(Probability pc, size_t idx) noexcept = 0;
        constexpr virtual Probability mutation_rate_impl(size_t idx) const noexcept = 0;

        constexpr virtual bool allow_variable_length_impl(size_t idx) const noexcept = 0;

        virtual void* component_impl(size_t idx) noexcept = 0;
        virtual const void* component_impl(size_t idx) const noexcept = 0;
    };

} // namespace gapp::mutation

#endif // !GAPP_MUTATION_BASE_DECL_HPP
