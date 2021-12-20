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


/** Implementations of commonly used crossover operators for integer encoded genetic algorithms. */

#ifndef GA_CROSSOVER_INTEGER_HPP
#define GA_CROSSOVER_INTEGER_HPP

#include "crossover_base.decl.hpp"
#include "../candidate.h"

#include <cstddef>


/** Predefined crossover operators for the integer encoded genetic algorithms (IntegerGA). */
namespace genetic_algorithm::crossover::integer
{
    /**
    * Standard single-point crossover operator for the integer encoded algorithms.
    * A random crossover point (locus) is selected and the genes before the locus are swapped
    * between the parents to create the children.
    */
    class SinglePoint final : public Crossover<size_t>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<size_t> crossover(const GA<size_t>& ga, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
    };

    /**
    * Two-point crossover operator for the integer encoded algorithms.
    * Two random crossover points are selected, and the genes between the 2 point are
    * swapped between the parents in order to create the children.
    * (Works as if 2 consecutive single-point crossovers were performed on the parents.)
    */
    class TwoPoint final : public Crossover<size_t>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<size_t> crossover(const GA<size_t>& ga, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
    };

    /**
    * Uniform crossover operator for the integer encoded algorithms.
    * Every gene of the chromosomes is swapped between the parents with 0.5 probability when creating the children.
    */
    class Uniform final : public Crossover<size_t>
    {
    public:
        using Crossover::Crossover;
    private:
        CandidatePair<size_t> crossover(const GA<size_t>& ga, const Candidate<size_t>& parent1, const Candidate<size_t>& parent2) const override;
    };

} // namespace genetic_algorithm::crossover::integer

#endif // !GA_CROSSOVER_INTEGER_HPP