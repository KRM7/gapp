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

/** Implementations of some commonly used mutation operators for integer encoded genetic algorithms. */

#ifndef GA_MUTATION_INTEGER_HPP
#define GA_MUTATION_INTEGER_HPP

#include "mutation_base.hpp"

#include <cstddef>

/** Predefined mutation operators for the integer encoded genetic algorithm. */
namespace genetic_algorithm::mutation::integer
{
    /**
    * Uniform mutation operator for the integer encoded genetic algorithm.
    * Each gene of the chromosome is changed with the specified mutation probability
    * to another value selected from a uniform distribution over all other values.
    */
    class Uniform : public Mutation<size_t>
    {
    public:
        /**
        * Create a uniform crossover operator with the specified parameters.
        * 
        * @param base The number of values a gene can take. Genes can be integers in the range [0, base-1].
        * @param The mutation probability used.
        */
        Uniform(size_t base, double pm = 0.01);

        /**
        * Sets the base for this crossover operator.
        * 
        * @param base The number of values a gene can take. Genes can be integers in the range [0, base-1].
        */
        void base(size_t base);

        /** @returns The base parameter currently set for this operator. */
        [[nodiscard]]
        size_t base() const noexcept { return base_; }

    private:
        void mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const override;

        size_t base_;
    };

} // namespace genetic_algorithm::mutation::integer

#endif // !GA_MUTATION_INTEGER_HPP