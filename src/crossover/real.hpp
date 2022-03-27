/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_REAL_HPP
#define GA_CROSSOVER_REAL_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"

#include <vector>
#include <utility>

/** Predefined crossover operators for the real encoded genetic algorithms (RCGA). */
namespace genetic_algorithm::crossover::real
{
    /**
    * Arithmetic crossover operator for the RCGA.
    * The children are the linear combinations of the parents, such that:
    *   child1 = alpha * parent1 + (1 - alpha) * parent2
    *   child2 = (1 - alpha) * parent1 + alpha * parent2
    * where alpha is a random number generated from a uniform distribution on [0.0, 1.0).
    */
    class Arithmetic final : public BoundedCrossover<double>
    {
    public:
        using BoundedCrossover::BoundedCrossover;
    private:
        CandidatePair<double> crossover(const GA<double>& ga, const Candidate<double>& parent1, const Candidate<double>& parent2) const override;
    };

    /**
    * BLX-Alpha (blend) crossover operator for the RCGA.
    * The genes of the children are chosen randomly from a uniform distribution based on
    * the values of the same genes of the parents.
    * 
    * The intervals the child genes are chosen from are:
    *   [-alpha * I + min(p1, p2), max(p1, p2) + alpha * I]
    * where I = abs(p1-p2)
    * 
    * This crossover operator has 1 parameter (alpha), which controls the length of the
    * intervals the child genes are chosen from. Larger alpha values -> larger intervals.
    * The recommended value of alpha is around 0.5.
    */
    class BLXa final : public BoundedCrossover<double>
    {
    public:
        /**
        * Create a BLX-alpha crossover operator with the specified parameters.
        *
        * @param bounds The (lower and upper) bounds of each gene. Must be the same length as the chromosomes used in the algorithm.
        * @param pc The crossover probability used.
        * @param alpha The alpha parameter of the BLX-alpha crossover. Must be >= 0.0.
        */
        explicit BLXa(const std::vector<std::pair<double, double>>& bounds, double pc = 0.8, double alpha = 0.5);

        /**
        * Sets the alpha parameter for the crossover.
        *
        * @param alpha The alpha parameter of the BLX-alpha crossover. Must be >= 0.0.
        */
        void alpha(double alpha);

        /** @returns The alpha parameter currently set for this operator. */
        [[nodiscard]]
        double alpha() const noexcept { return alpha_; }

    private:
        CandidatePair<double> crossover(const GA<double>& ga, const Candidate<double>& parent1, const Candidate<double>& parent2) const override;

        double alpha_ = 0.5;        /* The range parameter (alpha) of the BLX-alpha crossover. */
    };

    /**
    * Simulated binary crossover (SBX) operator for the RCGA.
    * The operator is based on the single-point crossover used in the binary encoded algorithms.
    * 
    * This crossover operator has one parameter, eta, which controls the shape of the probability
    * distribution the child genes are picked from. Larger eta values lead to children further away
    * from the parents, while smaller values will result in the children being closer to the parents.
    * Typical values for eta are [1.0, 5.0].
    */
    class SimulatedBinary final : public BoundedCrossover<double>
    {
    public:
        /**
        * Create a simulated binary crossover operator with the specified parameters.
        *
        * @param bounds The (lower and upper) bounds of each gene. Must be the same length as the chromosomes used in the algorithm.
        * @param pc The crossover probability used.
        * @param eta The shape parameter of the simulated binary crossover.
        */
        explicit SimulatedBinary(const std::vector<std::pair<double, double>>& bounds, double pc = 0.8, double eta = 4.0);

        /**
        * Sets the shape parameter (eta) of the simulated binary crossover.
        * 
        * @param eta The eta parameter of the crossover. Must be >=0.0.
        */
        void eta(double eta);

        /** @returns The eta parameter currently set for this operator. */
        [[nodiscard]]
        double eta() const noexcept { return eta_; }

    private:
        CandidatePair<double> crossover(const GA<double>& ga, const Candidate<double>& parent1, const Candidate<double>& parent2) const override;

        double eta_ = 4.0;          /* The shape parameter (eta) for the distribution of the simulated binary crossover. */
    };

    /**
    * Wright's heuristic crossover operator for the RCGA.
    * If p1 is the better parent, then the created children are:
    *   child1 = p1 + w1 * (p1 - p2)
    *   child2 = p1 + w2 * (p1 - p2)
    * where w1 and w2 are random numbers generated from a uniform distribution on [0.0, 1.0).
    */
    class Wright final : public BoundedCrossover<double>
    {
    public:
        using BoundedCrossover::BoundedCrossover;
    private:
        CandidatePair<double> crossover(const GA<double>& ga, const Candidate<double>& parent1, const Candidate<double>& parent2) const override;
    };


} // namespace genetic_algorithm::crossover::real

#endif // !GA_CROSSOVER_REAL_HPP