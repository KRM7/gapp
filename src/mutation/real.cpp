/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../population/candidate.hpp"
#include "../core/ga_base.hpp"
#include "../utility/rng.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>
#include <cstddef>

namespace genetic_algorithm::mutation::real
{
    void Uniform::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        GA_ASSERT(ga.gene_bounds().size() == candidate.chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(candidate.chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            candidate.chromosome[idx] = rng::randomReal(bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void NonUniform::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        GA_ASSERT(ga.gene_bounds().size() == candidate.chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(candidate.chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            const GeneType rand = rng::randomReal<GeneType>();
            const GeneType exponent = std::pow(1.0 - GeneType(ga.generation_cntr()) / ga.max_gen(), beta_);

            const GeneType multiplier = 1.0 - std::pow(rand, exponent);
            const GeneType bound = rng::randomBool() ? bounds[idx].lower() : bounds[idx].upper();

            candidate.chromosome[idx] += (bound - candidate.chromosome[idx]) * multiplier;
            /* The value of the mutated gene might be outside of the allowed interval. */
            candidate.chromosome[idx] = std::clamp(candidate.chromosome[idx], bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void Gauss::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        GA_ASSERT(ga.gene_bounds().size() == candidate.chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(candidate.chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            const GeneType SD = (bounds[idx].upper() - bounds[idx].lower()) / sigma_;

            candidate.chromosome[idx] += rng::randomNormal<GeneType>(0.0, SD);
            /* The value of the mutated gene might be outside of the allowed interval. */
            candidate.chromosome[idx] = std::clamp(candidate.chromosome[idx], bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void Polynomial::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        GA_ASSERT(ga.gene_bounds().size() == candidate.chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(candidate.chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            const GeneType alpha = rng::randomReal<GeneType>();
            if (alpha <= 0.5)
            {
                const GeneType delta = std::pow(2.0 * alpha, -(1.0 + eta_)) - 1.0;
                candidate.chromosome[idx] += delta * (candidate.chromosome[idx] - bounds[idx].lower());
            }
            else
            {
                const GeneType delta = 1.0 - std::pow(2.0 - 2.0 * alpha, -(1.0 + eta_));
                candidate.chromosome[idx] += delta * (bounds[idx].upper() - candidate.chromosome[idx]);
            }
            /* The value of the mutated gene might be outside of the allowed interval. */
            candidate.chromosome[idx] = std::clamp(candidate.chromosome[idx], bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void Boundary::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        GA_ASSERT(ga.gene_bounds().size() == candidate.chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(candidate.chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            candidate.chromosome[idx] = rng::randomBool() ? bounds[idx].lower() : bounds[idx].upper();
        }
    }

} // namespace genetic_algorithm::mutation::real