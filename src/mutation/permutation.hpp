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

/** Implementations of some commonly used mutation operators for the permutation encoded genetic algorithms. */

#ifndef GA_MUTATION_PERMUTATION_HPP
#define GA_MUTATION_PERMUTATION_HPP

#include "mutation_base.hpp"

/** Predefined mutation operators for the permutation encoded genetic algorithm. */
namespace genetic_algorithm::mutation::perm
{
    /**
    * Inversion mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the specified mutation probability.
    * In the mutated individuals, a random range of genes is chosen and reversed.
    */
    class Inversion : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Single swap/swap2 mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation probability.
    * Two randomly selected genes are swapped in the mutated individuals.
    */
    class Swap2 final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Swap-3 mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation probability.
    * Three randomly selected genes are swapped in the mutated individuals (0-1-2)->(2-0-1).
    */
    class Swap3 final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Shuffle/scramble mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation probability.
    * In the mutated candidates, a random range of genes is selected and then shuffled.
    */
    class Shuffle final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

    /**
    * Shift/slide mutation operator for the permutation encoded genetic algorithm.
    * Each individual is mutated with the set mutation proability.
    * In the mutated candidates, a random range of genes is selected and then moved to a
    * different position in the chromosome.
    */
    class Shift final : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::perm

#endif // !GA_MUTATION_PERMUTATION_HPP