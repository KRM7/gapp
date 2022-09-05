/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../encoding/integer.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>

namespace genetic_algorithm::mutation::integer
{
    void Uniform::mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const
    {
        const size_t base = dynamic_cast<const IntegerGA&>(ga).base();
        const size_t chrom_len = candidate.chromosome.size();

        const size_t mutate_count = rng::randomBinomial(chrom_len, mutation_rate());
        const auto mutated_indices = rng::sampleUnique(0_sz, chrom_len, mutate_count);

        if (base < 16) /* Make sure the new value for the changed genes can't be the old ones. */
        {
            std::vector<GeneType> alleles(base);
            std::iota(alleles.begin(), alleles.end(), GeneType{ 0 });

            for (const auto& idx : mutated_indices)
            {
                using std::swap;
                const auto old_gene = candidate.chromosome[idx];
                swap(alleles[old_gene], alleles.back());
                const GeneType new_gene = rng::randomElement(alleles.begin(), alleles.end() - 1);
                swap(alleles[old_gene], alleles.back());

                candidate.chromosome[idx] = new_gene;
            }
        }
        else /* The mutated gene values might be the same as the old ones, but the probability of this is low. */
        {
            for (const auto& idx : mutated_indices)
            {
                candidate.chromosome[idx] = rng::randomInt<GeneType>(0, base - 1);
            }
        }
    }

} // namespace genetic_algorithm::mutation::integer