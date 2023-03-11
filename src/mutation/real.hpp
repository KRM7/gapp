/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_REAL_HPP
#define GA_MUTATION_REAL_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/bounded_value.hpp"

/** Predefined mutation operators for the real encoded genetic algorithm (RCGA). */
namespace genetic_algorithm::mutation::real
{
    /**
    * Uniform mutation operator for the real encoded genetic algorithm (RCGA). \n
    * Each gene of the candidate is mutated with pm probability, and the values
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
    * Michalewicz's non-uniform mutation operator for the real encoded genetic algorithm (RCGA). \n
    * Each gene of the candidate is mutated with pm probability, and the values
    * of the mutated genes are randomly generated from a non-uniform distribution that changes over
    * time. In the early generations the distribution is close to uniform, while in the later generations
    * the mutated values tend to be closer to the original values. \n
    * 
    * The operator has one parameter, beta, which controls how fast the shape of the probability distribution
    * and how fast it changes over the generations. The value of this parameter must be >= 0.0. \n
    * For smaller values the distribution is more uniform and changes less over time (for beta = 0,
    * the distribution is uniform and doesn't change), while larger values will lead to faster change and
    * mutated genes closer to the original ones.
    */
    class NonUniform final : public Mutation<RealGene>
    {
    public:
        /**
        * Create a non-uniform mutation operator with the specified parameters.
        * 
        * @param pm The mutation probability used. Must be in the closed range [0.0, 1.0].
        * @param beta The beta parameter of the non-uniform crossover. Must be a finite value greater than 0.0.
        */
        constexpr explicit NonUniform(Probability pm, NonNegative<GeneType> beta = 2.0) noexcept :
            Mutation(pm), beta_(beta)
        {}
        
        /**
        * Sets the beta parameter for the crossover.
        * 
        * @param beta The beta parameter of the non-uniform crossover. Must be a finite value greater than 0.0.
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
    * Gauss mutation operator for the real encoded genetic algorithm (RCGA). \n
    * Each gene of the candidate is mutated with pm probability, and the values of the mutated genes
    * are randomly generated from a normal distribution around the current value of the gene. \n
    * 
    * The operator has one parameter, sigma, which controls the standard deviation of the normal distribution
    * used (but isn't the actual standard deviation). The SD of the normal distribution used for a gene is: \n
    *   SD = (upper_bound - lower_bound) / sigma \n
    * Larger sigma values will lead to mutated gene values closer to their original values.
    */
    class Gauss final : public Mutation<RealGene>
    {
    public:
        /**
        * Create a Gauss mutation operator with the specified parameters.
        * 
        * @param pm The mutation probability used. Must be in the closed range [0.0, 1.0].
        * @param sigma The sigma parameter of the gauss crossover. Must be a finite value greater than 0.0.
        */
        constexpr explicit Gauss(Probability pm, NonNegative<GeneType> sigma = 6.0) noexcept :
            Mutation(pm), sigma_(sigma)
        {}

        /**
        * Sets the sigma parameter for the crossover.
        * 
        * @param sigma The sigma parameter of the gauss crossover. Must be a finite value greater than 0.0.
        */
        constexpr void sigma(NonNegative<GeneType> sigma) noexcept { sigma_ = sigma; }

        /** @returns The sigma parameter of the operator. */
        [[nodiscard]]
        constexpr GeneType sigma() const noexcept { return sigma_; }

    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;

        NonNegative<GeneType> sigma_;
    };

    /**
    * Polynomial mutation operator for the real encoded genetic algorithm (RCGA). \n
    * Each gene of the candidate is mutated with pm probability, and the values of the mutated genes
    * are randomly generated from a non-uniform distribution. \n
    * 
    * This operator has one parameter, eta, which controls shape of the probability distribution
    * the mutated genes are picked from. The value of eta must be >= 0.0, with larger values leading
    * to mutated genes closer to the original ones. Typical values for eta are [20.0, 100.0].
    */
    class Polynomial final : public Mutation<RealGene>
    {
    public:
        /**
        * Create a polynomial mutation operator with the specified parameters.
        * 
        * @param bounds The (lower and upper) bounds of each gene. Must be the same length as the chromosomes used in the algorithm.
        * @param pm The mutation probability used. Must be in the closed range [0.0, 1.0].
        * @param eta The eta parameter of the polynomial mutation. Must be a finite value greater than 0.0.
        */
        constexpr explicit Polynomial(Probability pm, NonNegative<GeneType> eta = 40.0) noexcept :
            Mutation(pm), eta_(eta)
        {}

        /**
        * Sets the eta parameter for the crossover.
        * 
        * @param eta The eta parameter of the polynomial mutation. Must be a finite value greater than 0.0.
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
    * Boundary mutation operator for the real encoded genetic algorithm. \n
    * Each gene of the candidate is mutated with pm probability, and the values of the
    * mutated genes are either the lower or upper bounds of the gene.
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