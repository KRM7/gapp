/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BASE_DECL_HPP
#define GA_CROSSOVER_BASE_DECL_HPP

#include "../population/candidate.hpp"
#include <vector>
#include <utility>

namespace genetic_algorithm
{
    template<Gene T>
    class GA;
}

/** Crossover operators used in the algorithms. */
namespace genetic_algorithm::crossover
{
    /**
    * Base class used for all of the crossover operators.
    * Every implemented crossover operator takes 2 Candidates (the parents), and creates 2 children from these parents.
    * The crossover operation is performed on the 2 parents with the set (pc) probability only, the rest of the time
    * the returned children will be the same as the parents.
    */
    template<Gene T>
    class Crossover
    {
    public:
        /**
        * Create a crossover operator that will use @p pc as the crossover rate.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        explicit Crossover(double pc = 0.8);

        virtual ~Crossover() = default;

        /**
        * Sets the crossover rate used in the algorithm to @p pc. 
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        void crossover_rate(double pc);

        /** @returns The crossover rate currently set for this crossover operator. */
        [[nodiscard]]
        double crossover_rate() const noexcept { return pc_; }

        /**
        * Perform the crossover operation on 2 individuals with the set probability.
        *
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param parent1 The first parent.
        * @param parent2 The second parent.
        * @returns The pair of children resulting from the crossover.
        */
        CandidatePair<T> operator()(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const;

    protected:

        /* The actual crossover function. Performs the crossover every time and does nothing else. */
        virtual CandidatePair<T> crossover(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const = 0;

        double pc_ = 0.8;   /* The probability of performing the crossover operation on the parents. */
    };

    /**
    * Base class used for the crossover operators where each gene has a lower and upper bound.
    * Used for the operators that perform crossover on real encoded chromosomes.
    */
    template<Gene T>
    class BoundedCrossover : public Crossover<T>
    {
    public:
        /**
        * Create a bounded crossover operator with @p bounds used as the bounds for the genes.
        *
        * @param bounds The (lower and upper) bounds of the genes. Must be the same length as the chromosomes used in the algorithm.
        * @param pc The crossover probability.
        */
        explicit BoundedCrossover(const std::vector<std::pair<T, T>>& bounds, double pc = 0.8);

        /**
        * Sets the boundaries of the genes.
        * Each element of the @p bounds vector must contain the lower and upper bounds of the corresponding gene
        * (the min and max values of the gene), with the lower bound never being greater than the upper bound for a gene.
        * The number of elements in @p bounds must be the same as the length of the chromosomes used.
        *
        * @param bounds The (lower and upper) bounds of the genes.
        */
        void bounds(const std::vector<std::pair<T, T>>& bounds);

        /** @returns The bounds currently set for this crossover operator. */
        [[nodiscard]]
        std::vector<std::pair<T, T>> bounds() const noexcept { return bounds_; }

    protected:
        std::vector<std::pair<T, T>> bounds_;     /* The lower and upper bounds for each gene of the chromosome. */
    };

} // namespace genetic_algorithm::crossover

#endif // !GA_CROSSOVER_BASE_DECL_HPP