/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../utility/rng.hpp"
#include <vector>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <random>

namespace genetic_algorithm::mutation::binary
{
    void Flip::mutate(const GaInfo&, Candidate<char>& candidate) const
    {
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
            candidate.chromosome[idx] = !bool(candidate.chromosome[idx]);
        }
    }

} // namespace genetic_algorithm::mutation::binary