/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_REAL_HPP
#define GA_MUTATION_REAL_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"

/** Predefined mutation operators for the real encoded genetic algorithm. */
namespace genetic_algorithm::mutation::real
{
    /**
    * %Uniform mutation operator for the real encoded %GA.
    * 
    * Each gene of a candidate solution is mutated with the specified probability, and the values
    * of the mutated genes are randomly generated from a uniform distribution within the gene bounds.
    */
    class Uniform final : public Mutation<RealGene>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;
    };

    /**
    * Michalewicz's non-uniform mutation operator for the real encoded %GA.
    * 
    * Each gene of a candidate solution is mutated with the specified probability,
    * and the values of the mutated genes are randomly selected from a non-uniform
    * distribution that changes over time.
    * In the early generations the distribution is close to being uniform, while in the
    * later generations the mutated values tend to be closer to the original values.
    * 
    * The operator has one parameter (@p beta), which controls how quickly the shape of the
    * probability distribution changes over the generations. For smaller values, the distribution
    * is more uniform and changes less over time, while for larger values it will change faster
    * and the mutated genes will tend to be closer to the original ones.
    * 
    * The value of this parameter must be >= 0.0. If the value is 0, the distribution is
    * uniform and won't change at all over the generations.
    */
    class NonUniform final : public Mutation<RealGene>
    {
    public:
        /**
        * Create a non-uniform mutation operator with the specified parameters.
        * 
        * @param pm The mutation probability used. Must be in the closed range [0.0, 1.0].
        * @param beta The beta parameter of the mutation. Must be a finite non-negative value.
        */
        constexpr explicit NonUniform(Probability pm, NonNegative<GeneType> beta = 2.0) noexcept :
            Mutation(pm), beta_(beta)
        {}
        
        /**
        * Set the beta parameter of the mutation.
        * 
        * The value of this parameter controls how quickly the shape of the
        * probability distribution changes over the generations.
        * For smaller values, the distribution is more uniform and changes less over time,
        * while for larger values it will change faster and the mutated genes will
        * tend to be closer to the original ones.
        * 
        * @param beta The beta parameter of the mutation. Must be a finite non-negative value.
        */
        constexpr void beta(NonNegative<GeneType> beta) noexcept { beta_ = beta; }

        /** @returns The beta parameter of the operator. */
        [[nodiscard]]
        constexpr GeneType beta() const noexcept { return beta_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        NonNegative<GeneType> beta_;
    };

    /**
    * %Gauss mutation operator for the real encoded %GA.
    * 
    * Each gene of a candidate solution is mutated with the specified probability,
    * and the values of the mutated genes are randomly generated from a normal distribution
    * around the current value of the gene.
    * 
    * The operator has one parameter (@p sigma), which controls the standard deviation of the
    * normal distribution used (but isn't the actual standard deviation).
    * The standard deviation of the normal distribution used for the i-th gene is calculated as:
    *   \f[ SD_i = \frac{ bound_i^{upper} - bound_i^{lower} }{ sigma } \f]
    * 
    * Larger sigma values will lead to the values of the mutated gene to be closer to their
    * original values.
    */
    class Gauss final : public Mutation<RealGene>
    {
    public:
        /**
        * Create a %Gauss mutation operator with the specified parameters.
        * 
        * @param pm The mutation probability used. Must be in the closed range [0.0, 1.0].
        * @param sigma The sigma parameter of the crossover. Must be a finite value greater than 0.0.
        */
        constexpr explicit Gauss(Probability pm, Positive<GeneType> sigma = 6.0) noexcept :
            Mutation(pm), sigma_(sigma)
        {}

        /**
        * Set the sigma parameter of the mutation.
        * 
        * The value of this parameter controls the standard deviation of the normal distribution
        * the mutated gene values are selected from. The standard deviations are calculated from
        * the value of the parameter as:
        *   \f[ SD_i = \frac{ bound_i^{upper} - bound_i^{lower} }{ sigma } \f]
        * 
        * Larger sigma values will lead to the values of the mutated genes to be closer to their
        * original values.
        * 
        * @param sigma The sigma parameter of the gauss crossover. Must be a finite value greater than 0.0.
        */
        constexpr void sigma(Positive<GeneType> sigma) noexcept { sigma_ = sigma; }

        /** @returns The sigma parameter of the operator. */
        [[nodiscard]]
        constexpr GeneType sigma() const noexcept { return sigma_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        Positive<GeneType> sigma_;
    };

    /**
    * %Polynomial mutation operator for the real encoded %GA.
    * Each gene of a candidate solution is mutated with the specified probability,
    * and the values of the mutated genes are randomly generated from a non-uniform distribution.
    * 
    * This mutation operator has one parameter (@p eta), which controls shape of the probability distribution
    * the mutated genes are picked from. Larger values will lead to the mutated genes being closer
    * to their original values, while smaller ones will have the opposite effect.
    * 
    * The value of the parameter must be a finite, non-negative value.
    * Typical values are in the range [20.0, 100.0].
    */
    class Polynomial final : public Mutation<RealGene>
    {
    public:
        /**
        * Create a polynomial mutation operator with the specified parameters.
        * 
        * @param pm The mutation probability used. Must be in the closed range [0.0, 1.0].
        * @param eta The eta parameter of the mutation. Must be a finite non-negative value.
        */
        constexpr explicit Polynomial(Probability pm, NonNegative<GeneType> eta = 40.0) noexcept :
            Mutation(pm), eta_(eta)
        {}

        /**
        * Set the eta parameter for the crossover.
        * 
        * The value of this parameter controls the shape of the probability distribution
        * the mutated genes are picked from. Larger values will lead to the mutated genes
        * being closer to their original values, while using smaller values for the parameter
        * will have the opposite effect.
        * 
        * Typical values are in the range [20.0, 100.0].
        * 
        * @param eta The eta parameter of the mutation. Must be a finite non-negative value.
        */
        constexpr void eta(NonNegative<GeneType> eta) noexcept { eta_ = eta; }

        /** @returns The eta parameter of the operator. */
        [[nodiscard]]
        constexpr GeneType eta() const noexcept { return eta_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        NonNegative<GeneType> eta_;
    };

    /**
    * %Boundary mutation operator for the real encoded %GA.
    * 
    * Each gene of a candidate solution is mutated with the specified probability,
    * and the values of the mutated genes are either the lower or upper bounds of the gene
    * (with equal probability).
    */
    class Boundary final : public Mutation<RealGene>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::real

#endif // !GA_MUTATION_REAL_HPP