/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/bounded_value.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <tuple>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::mutation::perm
{
    void Inversion::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        const size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() < mutation_rate())
        {
            const size_t min_len = 2;
            const size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            const size_t range_len = rng::randomInt(min_len, max_len);

            const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
            const size_t last  = first + range_len;

            std::reverse(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last);
        }
    }

    void Swap2::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        if (candidate.chromosome.size() < 2) return;

        if (rng::randomReal() <= mutation_rate())
        {
            const auto idxs = rng::sampleUnique(0_sz, candidate.chromosome.size(), 2_sz);

            using std::swap;
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
        }
    }

    void Swap3::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        if (candidate.chromosome.size() < 3) return;

        if (rng::randomReal() <= mutation_rate())
        {
            const auto idxs = rng::sampleUnique(0_sz, candidate.chromosome.size(), 3_sz);

            using std::swap;
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[2]]);
        }
    }

    void Shuffle::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        const size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() <= mutation_rate())
        {
            const size_t min_len = 2;
            const size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            const size_t range_len = rng::randomInt(min_len, max_len);

            const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
            const size_t last  = first + range_len;

            std::shuffle(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last, rng::prng);
        }
    }

    void Shift::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        const size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() <= mutation_rate())
        {
            const size_t min_len = 1;
            const size_t max_len = std::max(min_len, size_t(range_max_ * chrom_len));
            const size_t range_len = rng::randomInt(min_len, max_len);

            const auto idxs = rng::sampleUnique(0_sz, chrom_len - range_len + 1, 2);
            const auto src_first  = candidate.chromosome.begin() + ptrdiff_t(idxs[0]);
            const auto dest_first = candidate.chromosome.begin() + ptrdiff_t(idxs[1]);

            const auto [first, middle, last] = (dest_first < src_first) ?
                std::tuple(dest_first, src_first, src_first + range_len) :
                std::tuple(src_first, src_first + range_len, dest_first + range_len);

            std::rotate(first, middle, last);
        }
    }

} // namespace genetic_algorithm::mutation::perm