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

#include "binary.hpp"
#include "../candidate.h"
#include "../utils.h"
#include "../rng.hpp"

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cassert>
#include <cstddef>

namespace genetic_algorithm::crossover::binary
{
    CandidatePair<char> nPointCrossoverImpl(const Candidate<char>& parent1, const Candidate<char>& parent2, size_t n)
    {
        assert(n > 0);

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the n-point crossover.");
        }

        std::vector<size_t> loci = rng::sampleUnique(parent1.chromosome.size(), std::min(n, parent1.chromosome.size()));

        /* Count how many loci are after each gene. */
        std::vector<size_t> loci_after;
        loci_after.reserve(parent1.chromosome.size());

        for (size_t i = 0, loci_left = loci.size(); i < parent1.chromosome.size(); i++)
        {
            if (loci_left > 0 && std::find(loci.begin(), loci.end(), i) != loci.end())
            {
                loci_left--;
            }
            loci_after.push_back(loci_left);
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            if (loci_after[i] % 2)
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }

        return { child1, child2 };
    }

    CandidatePair<char> SinglePoint::crossover(const GA<char>& ga, const Candidate<char>& parent1, const Candidate<char>& parent2) const
    {
        GA_UNUSED(ga);

        return nPointCrossoverImpl(parent1, parent2, 1U);
    }

    CandidatePair<char> TwoPoint::crossover(const GA<char>& ga, const Candidate<char>& parent1, const Candidate<char>& parent2) const
    {
        GA_UNUSED(ga);

        return nPointCrossoverImpl(parent1, parent2, 2U);
    }

    CandidatePair<char> Uniform::crossover(const GA<char>& ga, const Candidate<char>& parent1, const Candidate<char>& parent2) const
    {
        GA_UNUSED(ga);

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the uniform crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            /* Swap each gene with 0.5 probability. */
            if (rng::randomBool())
            {
                child1.chromosome[i] = parent2.chromosome[i];
                child2.chromosome[i] = parent1.chromosome[i];
            }
        }

        return { child1, child2 };
    }

} // namespace genetic_algorithm::crossover::binary