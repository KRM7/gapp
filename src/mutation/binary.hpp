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

/** Implementations of commonly used mutation operators for the binary encoded genetic algorithms. */

#ifndef GA_MUTATION_BINARY_HPP
#define GA_MUTATION_BINARY_HPP

#include "mutation_base.hpp"

namespace genetic_algorithm::mutation::binary
{
    /**
    * Standard flip mutation for the binary encoded genetic algorithm.
    * Each gene of the chromosome is flipped with the specified mutation probability.
    */
    class Flip final : public Mutation<char>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<char>& ga, Candidate<char>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::binary

#endif // !GA_MUTATION_BINARY_HPP