/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_REAL_HPP
#define GA_MUTATION_REAL_HPP

#include "mutation_base.hpp"

#include <vector>
#include <utility>

/** Predefined mutation operators for the real encoded genetic algorithm (RCGA). */
namespace genetic_algorithm::mutation::real
{
    /**
    * Uniform mutation operator for the real encoded genetic algorithm (RCGA).
    * Each gene of the candidate is mutated with pm probability, and the values
    * of the mutated genes are randomly generated from a uniform distribution within the gene bounds.
    */
    class Uniform final : public BoundedMutation<double>
    {
    public:
        using BoundedMutation::BoundedMutation;
    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

    /**
    * Michalewicz's non-uniform mutation operator for the real encoded genetic algorithm (RCGA).
    * Each gene of the candidate is mutated with pm probability, and the values
    * of the mutated genes are randomly generated from a non-uniform distribution that changes over
    * time. In the early generations the distribution is close to uniform, while in the later generations
    * the mutated values tend to be closer to the original values.
    * 
    * The operator has one parameter, beta, which controls how fast the shape of the probability distribution
    * and how fast it changes over the generations. The value of this parameter must be >= 0.0.
    * For smaller values the distribution is more uniform and changes less over time (for beta = 0,
    * the distribution is uniform and doesn't change), while larger values will lead to faster change and
    * mutated genes closer to the original ones.
    */
    class NonUniform final : public BoundedMutation<double>
    {
    public:
        /**
        * Create a non-uniform mutation operator with the specified parameters.
        * 
        * @param bounds The (lower and upper) bounds of each gene. Must be the same length as the chromosomes used in the algorithm.
        * @param pm The mutation probability used.
        * @param beta The beta parameter of the non-uniform crossover. Must be >= 0.0.
        */
        explicit NonUniform(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double beta = 2.0);
        
        /**
        * Sets the beta parameter for the crossover.
        * 
        * @param beta The beta parameter of the non-uniform crossover. Must be >= 0.0.
        */
        void beta(double beta);

        /** @returns The beta parameter currently set for this operator. */
        [[nodiscard]]
        double beta() const noexcept { return beta_; }

    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;

        double beta_ = 2.0;
    };

    /**
    * Gauss mutation operator for the real encoded genetic algorithm (RCGA).
    * Each gene of the candidate is mutated with pm probability, and the values of the mutated genes
    * are randomly generated from a normal distribution around the current value of the gene.
    * 
    * The operator has one parameter, sigma, which controls the standard deviation of the normal distribution
    * used (but isn't the actual standard deviation). The SD of the normal distribution used for a gene is:
    *   SD = (upper_bound - lower_bound) / sigma
    * Larger sigma values will lead to mutated gene values closer to their original values.
    */
    class Gauss final : public BoundedMutation<double>
    {
    public:
        /**
        * Create a Gauss mutation operator with the specified parameters.
        * 
        * @param bounds The (lower and upper) bounds of each gene. Must be the same length as the chromosomes used in the algorithm.
        * @param pm The mutation probability used.
        * @param sigma The sigma parameter of the gauss crossover. Must be > 0.0.
        */
        explicit Gauss(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double sigma = 6.0);

        /**
        * Sets the sigma parameter for the crossover.
        * 
        * @param sigma The sigma parameter of the gauss crossover. Must be > 0.0.
        */
        void sigma(double sigma);

        /** @returns The sigma parameter currently set for this operator. */
        [[nodiscard]]
        double sigma() const noexcept { return sigma_; }

    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;

        double sigma_ = 6.0;
    };

    /**
    * Polynomial mutation operator for the real encoded genetic algorithm (RCGA).
    * Each gene of the candidate is mutated with pm probability, and the values of the mutated genes
    * are randomly generated from a non-uniform distribution.
    * 
    * This operator has one parameter, eta, which controls shape of the probability distribution
    * the mutated genes are picked from. The value of eta must be >= 0.0, with larger values leading
    * to mutated genes closer to the original ones. Typical values for eta are [20.0, 100.0].
    */
    class Polynomial final : public BoundedMutation<double>
    {
    public:
        /**
        * Create a polynomial mutation operator with the specified parameters.
        * 
        * @param bounds The (lower and upper) bounds of each gene. Must be the same length as the chromosomes used in the algorithm.
        * @param pm The mutation probability used.
        * @param eta The eta parameter of the polynomial mutation. Must be >= 0.0.
        */
        explicit Polynomial(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double eta = 40.0);

        /**
        * Sets the eta parameter for the crossover.
        * 
        * @param eta The eta parameter of the polynomial mutation. Must be >= 0.0.
        */
        void eta(double eta);

        /** @returns The eta parameter currently set for this operator. */
        [[nodiscard]] double eta() const noexcept { return eta_; }

    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;

        double eta_ = 40.0;
    };

    /**
    * Boundary mutation operator for the real encoded genetic algorithm.
    * Each gene of the candidate is mutated with pm probability, and the values of the
    * mutated genes are either the lower or upper bounds of the given gene.
    */
    class Boundary final : public BoundedMutation<double>
    {
    public:
        using BoundedMutation::BoundedMutation;
    private:
        void mutate(const GA<double>& ga, Candidate<double>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::real

#endif // !GA_MUTATION_REAL_HPP