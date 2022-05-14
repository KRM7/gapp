/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "mutation_dtl.hpp"
#include "../algorithms/integer_ga.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <vector>

namespace genetic_algorithm::mutation::integer
{
    void Uniform::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        size_t base = dynamic_cast<const IntegerGA&>(ga).base();
        size_t chrom_len = candidate.chromosome.size();

        size_t mutate_count = dtl::approxMutateCnt(chrom_len, pm_);
        auto mutated_indices = rng::sampleUnique(0_sz, chrom_len, mutate_count);

        std::vector<GeneType> alleles(base);
        std::iota(alleles.begin(), alleles.end(), GeneType{ 0 });

        for (const auto& idx : mutated_indices)
        {
            /* Make sure the new value for the changed gene can't be the old one. */
            std::swap(alleles[candidate.chromosome[idx]], alleles.back());
            GeneType new_gene = rng::randomElement(alleles.begin(), std::prev(alleles.end()));
            std::swap(alleles[candidate.chromosome[idx]], alleles.back());

            candidate.chromosome[idx] = new_gene;
        }
    }

} // namespace genetic_algorithm::mutation::integer