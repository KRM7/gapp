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
#include "../utils.hpp"
#include "../rng.hpp"

#include <vector>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <random>

namespace genetic_algorithm::mutation::binary
{
    void Flip::mutate(const GA<char>& ga, Candidate<char>& candidate) const
    {
        GA_UNUSED(ga);

        size_t flip_cnt;
        if (candidate.chromosome.size() * pm_ >= 2.0)
        {
            double mean = candidate.chromosome.size() * pm_;
            double SD = std::sqrt(mean * (1.0 - pm_));

            double res = rng::randomNormal(mean, SD);
            while (res <= -0.5) { res = rng::randomNormal(mean, SD); }

            flip_cnt = std::min(size_t(std::round(res)), candidate.chromosome.size());
        }
        else
        {
            std::binomial_distribution<size_t> dist(candidate.chromosome.size(), pm_);
            flip_cnt = dist(rng::prng);
        }

        std::vector<size_t> flipped_idxs = rng::sampleUnique(candidate.chromosome.size(), flip_cnt);
        for (const auto& idx : flipped_idxs)
        {
            candidate.chromosome[idx] = !candidate.chromosome[idx];
        }
    }

} // namespace genetic_algorithm::mutation::binary