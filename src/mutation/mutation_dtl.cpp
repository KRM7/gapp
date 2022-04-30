/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "mutation_dtl.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <random>
#include <cmath>

namespace genetic_algorithm::mutation::dtl
{
    size_t approxMutateCnt(size_t chrom_len, double pm)
    {
        size_t mutate_cnt;
        if (chrom_len * pm >= 2.0)
        {
            double mean = chrom_len * pm;
            double SD = std::sqrt(mean * (1.0 - pm));

            double r = rng::randomNormal(mean, SD);
            while (r <= -0.5) { r = rng::randomNormal(mean, SD); }

            mutate_cnt = size_t(std::round(r));
            mutate_cnt = std::min(mutate_cnt, chrom_len);
        }
        else
        {
            std::binomial_distribution dist(chrom_len, pm);
            mutate_cnt = dist(rng::prng);
        }

        return mutate_cnt;
    }

} // namespace genetic_algorithm::mutation::dtl