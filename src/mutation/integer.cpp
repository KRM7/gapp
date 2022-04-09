/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../utility/rng.hpp"
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

    void Uniform::mutate(const GaInfo&, Candidate<size_t>& candidate) const
    {
        size_t mutate_cnt;
        if (candidate.chromosome.size() * pm_ >= 2.0)
        {
            double mean = candidate.chromosome.size() * pm_;
            double SD = std::sqrt(mean * (1.0 - pm_));

            double res = rng::randomNormal(mean, SD);
            while (res <= -0.5) { res = rng::randomNormal(mean, SD); }

            mutate_cnt = std::min(size_t(std::round(res)), candidate.chromosome.size());
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
            size_t new_gene = alleles[rng::randomInt(size_t{ 0 }, alleles.size() - 2)];
            std::swap(alleles[candidate.chromosome[idx]], alleles.back());

            candidate.chromosome[idx] = new_gene;
        }
    }

} // namespace genetic_algorithm::mutation::integer