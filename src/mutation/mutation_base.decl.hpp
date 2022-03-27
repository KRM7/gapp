/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_DECL_HPP
#define GA_MUTATION_BASE_DECL_HPP

#include "../candidate.hpp"
#include "../concepts.hpp"

#include <vector>
#include <utility>

namespace genetic_algorithm
{
    template<typename T>
    class GA;

} // namespace genetic_algorithm

/** Mutation operators used in the algorithms. */
namespace genetic_algorithm::mutation
{
    /**
    * Base class used for all of the mutation operators.
    * Every implemented mutation operator takes a Candidate, and mutates that candidate using the set probability.
    * (This probability can either be per-candidate or per-gene depending on the operator.)
    */
    template<Gene T>
    class Mutation
    {
    public:
        /**
        * Create a mutation operator that will use @p pm as the mutation rate.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        explicit Mutation(double pm = 0.01);

        virtual ~Mutation() = default;

        /**
        * Sets the mutation rate used for the operator to @p pm.
        * 
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        void mutation_rate(double pm);

        /** @returns The mutation rate currently set for this mutation operator. */
        [[nodiscard]]
        double mutation_rate() const noexcept { return pm_; }

        /**
        * Perform the mutation operation on a candidate using the set mutation rate.
        *
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param candidate The candidate to mutate.
        */
        void operator()(const GA<T>& ga, Candidate<T>& candidate) const;

    protected:

        /* The actual mutation function. Performs the mutation using the set probability and does nothing else. */
        virtual void mutate(const GA<T>& ga, Candidate<T>& candidate) const = 0;

        double pm_ = 0.01;      /* The mutation rate used for the mutation operator. */
    };

    /**
    * Base class used for the mutation operators where each gene also has a lower and upper bound.
    * Used for the operators that perform mutation on real encoded chromosomes.
    */
    template<Gene T>
    class BoundedMutation : public Mutation<T>
    {
    public:
        /**
        * Create a bounded mutation operator with @p bounds used as the bounds for the genes.
        *
        * @param bounds The (lower and upper) bounds of the genes. Must be the same length as the chromosomes used in the algorithm.
        * @param pm The mutation probability.
        */
        explicit BoundedMutation(const std::vector<std::pair<T, T>>& bounds, double pm = 0.01);

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

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_DECL_HPP