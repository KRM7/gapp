/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../core/candidate.hpp"
#include "../core/ga_base.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstddef>

namespace gapp::mutation::real
{
    void Uniform::mutate(const GA<GeneType>& ga, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        GAPP_ASSERT(ga.gene_bounds().size() == chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            chromosome[idx] = rng::randomReal(bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void NonUniform::mutate(const GA<GeneType>& ga, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        GAPP_ASSERT(ga.gene_bounds().size() == chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            const GeneType rand = rng::randomReal<GeneType>();
            const GeneType exponent = std::pow(1.0 - GeneType(ga.generation_cntr()) / ga.max_gen(), beta_);

            const GeneType multiplier = 1.0 - std::pow(rand, exponent);
            const GeneType bound = rng::randomBool() ? bounds[idx].lower() : bounds[idx].upper();

            chromosome[idx] += (bound - chromosome[idx]) * multiplier;
            /* The value of the mutated gene might be outside of the allowed interval. */
            chromosome[idx] = std::clamp(chromosome[idx], bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void Gauss::mutate(const GA<GeneType>& ga, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        GAPP_ASSERT(ga.gene_bounds().size() == chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            const GeneType SD = (bounds[idx].upper() - bounds[idx].lower()) / sigma_;

            chromosome[idx] += rng::randomNormal<GeneType>(0.0, SD);
            /* The value of the mutated gene might be outside of the allowed interval. */
            chromosome[idx] = std::clamp(chromosome[idx], bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void Polynomial::mutate(const GA<GeneType>& ga, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        GAPP_ASSERT(ga.gene_bounds().size() == chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            const GeneType alpha = rng::randomReal<GeneType>();
            if (alpha <= 0.5)
            {
                const GeneType delta = std::pow(2.0 * alpha, 1.0 / (1.0 + eta_)) - 1.0;
                chromosome[idx] += delta * (chromosome[idx] - bounds[idx].lower());
            }
            else
            {
                const GeneType delta = 1.0 - std::pow(2.0 - 2.0 * alpha, 1.0 / (1.0 + eta_));
                chromosome[idx] += delta * (bounds[idx].upper() - chromosome[idx]);
            }
            /* The value of the mutated gene might be outside of the allowed interval. */
            chromosome[idx] = std::clamp(chromosome[idx], bounds[idx].lower(), bounds[idx].upper());
        }
    }

    void Boundary::mutate(const GA<GeneType>& ga, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const
    {
        GAPP_ASSERT(ga.gene_bounds().size() == chromosome.size(), "Mismatching bounds and chromosome lengths.");

        const auto& bounds = ga.gene_bounds();

        const size_t mutate_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            chromosome[idx] = rng::randomBool() ? bounds[idx].lower() : bounds[idx].upper();
        }
    }

} // namespace gapp::mutation::real