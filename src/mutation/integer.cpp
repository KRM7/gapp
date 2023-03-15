/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../population/candidate.hpp"
#include "../core/ga_base.hpp"
#include "../utility/rng.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::mutation::integer
{
    void Uniform::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        const auto& bounds = ga.gene_bounds();
        const size_t chrom_len = candidate.chromosome.size();

        const size_t mutate_count = rng::randomBinomial(chrom_len, mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chrom_len, mutate_count);

        for (const auto& idx : mutated_indices)
        {
            GeneType old_gene = candidate.chromosome[idx];
            GeneType new_gene = rng::randomInt(bounds[idx].lower(), bounds[idx].upper());

            while (new_gene == old_gene)
            {
                new_gene = rng::randomInt(bounds[idx].lower(), bounds[idx].upper());
            }

            candidate.chromosome[idx] = new_gene;
        }
    }

} // namespace genetic_algorithm::mutation::integer