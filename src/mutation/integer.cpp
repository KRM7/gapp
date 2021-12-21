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

#include "integer.hpp"
#include "../utils.h"
#include "../rng.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

namespace genetic_algorithm::mutation::integer
{
    Uniform::Uniform(size_t base, double pm) :
        Mutation(pm)
    {
        this->base(base);
    }

    void Uniform::base(size_t base)
    {
        if (base < 2) throw std::invalid_argument("The base for the uniform integer mutation must be at least 2.");

        base_ = base;
    }

    void Uniform::mutate(const GA<size_t>& ga, Candidate<size_t>& candidate) const
    {
        GA_UNUSED(ga);

        size_t mutate_cnt;
        if (candidate.chromosome.size() * pm_ >= 2.0)
        {
            double mean = candidate.chromosome.size() * pm_;
            double SD = std::sqrt(mean * (1.0 - pm_));

            double res = rng::randomNormal(mean, SD);
            while (res <= -0.5) { res = rng::randomNormal(mean, SD); }

            mutate_cnt = std::min(size_t{ std::round(res) }, candidate.chromosome.size());
        }
        else
        {
            std::binomial_distribution<size_t> dist(candidate.chromosome.size(), pm_);
            mutate_cnt = dist(rng::prng);
        }

        std::vector<size_t> changed_idxs = rng::sampleUnique(candidate.chromosome.size(), mutate_cnt);

        std::vector<size_t> alleles(base_);
        std::iota(alleles.begin(), alleles.end(), size_t{ 0 });

        for (const auto& idx : changed_idxs)
        {
            /* Make sure the new value for the changed gene can't be the old one. */
            std::swap(alleles[candidate.chromosome[idx]], alleles.back());
            size_t new_gene = alleles[rng::randomIdx(alleles.size() - 1)];
            std::swap(alleles[candidate.chromosome[idx]], alleles.back());

            candidate.chromosome[idx] = new_gene;
        }
    }

} // namespace genetic_algorithm::mutation::integer