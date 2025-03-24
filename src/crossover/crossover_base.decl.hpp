/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_BASE_DECL_HPP
#define GAPP_CROSSOVER_BASE_DECL_HPP

#include "../core/candidate.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"
#include <array>
#include <utility>
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

namespace gapp::crossover
{
    /**
    * The base class used for the crossover operators of the GAs.
    * 
    * %Crossover operators take 2 candidate solutions (the parents), and create 2 new
    * candidates (children) based on the the parent candidates.
    * The crossover operation is only performed on the 2 parents with a set probability only,
    * the rest of the time the returned children will be the same as the parents.
    * 
    * New crossover operators should be derived from this class, and they must implement the
    * following virtual method:
    * 
    *   - crossover : Perform the crossover on 2 candidate solutions.
    * 
    * @tparam The gene type the crossover operator is defined for.
    */
    template<typename T>
    class Crossover
    {
    public:
        /** The gene type the crossover operator is defined for. */
        using GeneType = T;

        /** Create a crossover operator using the default crossover probability. */
        constexpr explicit Crossover() noexcept :
            pc_(0.8_p) {}

        /**
        * Create a crossover operator with the specified crossover probability.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Crossover(Probability pc) noexcept :
            pc_(pc) {}

        /**
        * Set the crossover probability used for the crossovers.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void crossover_rate(Probability pc) noexcept { pc_ = pc; }

        /** @returns The crossover rate set for the operator. */
        [[nodiscard]]
        constexpr Probability crossover_rate() const noexcept { return pc_; }

        /**
        * This method specifies whether the crossover operator supports variable
        * chromosome lengths or not. If variable chromosome lengths are supported,
        * the Candidates passed to the crossover operator are allowed to have
        * chromosome lengths different from each other, and from the chromosome length
        * specified for the GA the operator is used in. Otherwise the chromosome length
        * of every Candidate must be the same.
        *
        * This method will return false by default. If a particular crossover method allows
        * variable chromosome lengths, it should override this method to return true.
        *
        * @returns True if the crossover operator supports variable chromosome lengths.
        */
        [[nodiscard]]
        constexpr virtual bool allow_variable_chrom_length() const noexcept { return false; }

        /**
        * Perform the crossover operation on 2 candidate solutions with the set probability.
        * This function is implemented by crossover().
        *
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param parent1 The first parent solution.
        * @param parent2 The second parent solution.
        * @returns The pair of children resulting from the crossover.
        */
        CandidatePair<T> operator()(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const;


        /** Destructor. */
        virtual ~Crossover()                    = default;

    protected:

        Crossover(const Crossover&)             = default;
        Crossover(Crossover&&)                  = default;
        Crossover& operator=(const Crossover&)  = default;
        Crossover& operator=(Crossover&&)       = default;

    private:

        /**
        * The implementation of the crossover operator. Performs the crossover operation
        * on 2 parent solutions to generate 2 child solutions from them.
        * The implementation of this function shouldn't handle the crossover probability,
        * instead it should just perform the crossover operation unconditionally.
        * The chromosomes of the returned children should be valid solutions for the given
        * problem and %GA, but the rest of their properties (eg. fitness) are irrelevant.
        * 
        * This method will be called once for every 2 children that need to be generated
        * (ie. population_size/2 number of times, rounded up if the population size is odd)
        * in every generation.
        * 
        * The implementation of this function must be thread-safe.
        * 
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param parent1 The first parent solution.
        * @param parent2 The second parent solution.
        * @returns The pair of children resulting from the crossover.
        */
        virtual CandidatePair<T> crossover(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const = 0;

        Probability pc_;
    };

    /**
    * The base class used for the mixed crossover operators.
    * This is effectively the same as the primary template without the crossover
    * probability.
    *
    * @tparam Ts The component gene types of the mixed gene type the crossover operator is defined for.
    */
    template<typename... Ts>
    class Crossover<MixedGene<Ts...>>
    {
    public:
        /** The mixed gene type the crossover operator is defined for. */
        using GeneType = MixedGene<Ts...>;

        /** The number of component crossovers the mixed crossover is composed of. */
        static constexpr size_t N = sizeof...(Ts);

        /**
        * Set the crossover probability of the component crossover associated with the
        * specified GeneType.
        *
        * @tparam GeneType The gene type of the component crossover to set the probability of.
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        template<typename GeneType>
        constexpr void crossover_rate(Probability pc) noexcept;

        /**
        * Set the crossover probability used for each of the component crossovers to the
        * same value.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr virtual void crossover_rates(Probability pc) noexcept = 0;

        /**
        * Set the crossover probability used for each of the component crossovers
        * individually. The order of the probabilities should match the order of the
        * component crossovers.
        *
        * @param pcs The crossover probabilities. They must all be in the closed interval [0.0, 1.0].
        */
        constexpr virtual void crossover_rates(const std::array<Probability, N>& pcs) noexcept = 0;

        /**
        * @returns The crossover probability of the component crossover associated with the specified GeneType.
        * @tparam GeneType The gene type of the component crossover to get the probability of.
        */
        template<typename GeneType>
        [[nodiscard]] constexpr Probability crossover_rate() const noexcept;

        /**
        * @returns The crossover rates set for the component crossovers. The order of the
        *   probabilities in the returned array will match the order of the component crossovers.
        */
        [[nodiscard]]
        constexpr virtual std::array<Probability, N> crossover_rates() const noexcept = 0;

        /**
        * This method specifies whether the crossover operator supports variable
        * chromosome lengths for the specified gene type or not. If variable chromosome
        * lengths are supported, the Candidates passed to the crossover operator are
        * allowed to have chromosome lengths different from each other, and from the
        * chromosome length specified for the GA the operator is used in. Otherwise the
        * chromosome length of the given gene type must be the same for every Candidate.
        *
        * @returns True if the crossover operator supports variable chromosome lengths
        *   for the specified gene type.
        * @tparam GeneType The gene type of the chromosome to check.
        */
        template<typename GeneType>
        [[nodiscard]] constexpr bool allow_variable_chrom_length() const noexcept;

        /**
        * @returns The component crossover associated with the specified gene type.
        * @tparam GeneType The gene type of the component crossover to get.
        */
        template<typename GeneType>
        [[nodiscard]] Crossover<GeneType>& component() & noexcept;

        /**
        * @returns The component crossover associated with the specified gene type.
        * @tparam GeneType The gene type of the component crossover to get.
        */
        template<typename GeneType>
        [[nodiscard]] const Crossover<GeneType>& component() const& noexcept;

        /**
        * Perform the crossover operation on 2 candidate solutions.
        * This function is implemented by crossover().
        *
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param parent1 The first parent solution.
        * @param parent2 The second parent solution.
        * @returns The pair of children resulting from the crossover.
        */
        CandidatePair<GeneType> operator()(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const;


        /** Destructor. */
        virtual ~Crossover()                   = default;

    protected:

        Crossover()                            = default;
        Crossover(const Crossover&)            = default;
        Crossover(Crossover&&)                 = default;
        Crossover& operator=(const Crossover&) = default;
        Crossover& operator=(Crossover&&)      = default;

    private:
        virtual CandidatePair<GeneType> crossover(const GaInfo& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const = 0;

        constexpr virtual void crossover_rate_impl(Probability pc, size_t idx) noexcept = 0;
        constexpr virtual Probability crossover_rate_impl(size_t idx) const noexcept = 0;

        constexpr virtual bool allow_variable_length_impl(size_t idx) const noexcept = 0;

        virtual void* component_impl(size_t idx) noexcept = 0;
        virtual const void* component_impl(size_t idx) const noexcept = 0;
    };

} // namespace gapp::crossover

#endif // !GAPP_CROSSOVER_BASE_DECL_HPP
