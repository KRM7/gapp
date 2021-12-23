/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

#ifndef GA_MUTATION_BASE_DECL_HPP
#define GA_MUTATION_BASE_DECL_HPP

#include "../candidate.hpp"

#include <vector>
#include <utility>

namespace genetic_algorithm
{
    template<typename GeneType>
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
    template<regular_hashable GeneType>
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
        void operator()(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const;

    protected:

        /* The actual mutation function. Performs the mutation using the set probability and does nothing else. */
        virtual void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const = 0;

        double pm_ = 0.01;      /* The mutation rate used for the mutation operator. */
    };

    /**
    * Base class used for the mutation operators where each gene also has a lower and upper bound.
    * Used for the operators that perform mutation on real encoded chromosomes.
    */
    template<regular_hashable GeneType>
    class BoundedMutation : public Mutation<GeneType>
    {
    public:
        /**
        * Create a bounded mutation operator with @p bounds used as the bounds for the genes.
        *
        * @param bounds The (lower and upper) bounds of the genes. Must be the same length as the chromosomes used in the algorithm.
        * @param pm The mutation probability.
        */
        explicit BoundedMutation(const std::vector<std::pair<GeneType, GeneType>>& bounds, double pm = 0.01);

        /**
        * Sets the boundaries of the genes.
        * Each element of the @p bounds vector must contain the lower and upper bounds of the corresponding gene
        * (the min and max values of the gene), with the lower bound never being greater than the upper bound for a gene.
        * The number of elements in @p bounds must be the same as the length of the chromosomes used.
        *
        * @param bounds The (lower and upper) bounds of the genes.
        */
        void bounds(const std::vector<std::pair<GeneType, GeneType>>& bounds);

        /** @returns The bounds currently set for this crossover operator. */
        [[nodiscard]]
        std::vector<std::pair<GeneType, GeneType>> bounds() const noexcept { return bounds_; }

    protected:
        std::vector<std::pair<GeneType, GeneType>> bounds_;     /* The lower and upper bounds for each gene of the chromosome. */
    };

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_DECL_HPP