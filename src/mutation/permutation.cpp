/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/rng.hpp"
#include "../utility/probability.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cstddef>

namespace genetic_algorithm::mutation::perm
{
    Inversion::Inversion(Probability pm, double range_max)
        : Mutation(pm)
    {
        this->range_max(range_max);
    }

    void Inversion::range_max(double rm)
    {
        if (!(0.0 <= rm && rm <= 1.0))
        {
            GA_THROW(std::invalid_argument, "The range_max parameter must be in the closed interval [0.0, 1.0].");
        }

        range_max_ = rm;
    }

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
            const size_t last = first + range_len;

            std::reverse(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last);
        }
    }

    void Swap2::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        if (candidate.chromosome.size() < 2) return;

        if (rng::randomReal() < mutation_rate())
        {
            const auto idxs = rng::sampleUnique(0_sz, candidate.chromosome.size(), 2_sz);

            using std::swap;
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
        }
    }

    void Swap3::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        if (candidate.chromosome.size() < 3) return;

        if (rng::randomReal() < mutation_rate())
        {
            const auto idxs = rng::sampleUnique(0_sz, candidate.chromosome.size(), 3_sz);

            using std::swap;
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[2]]);
        }
    }

    Shuffle::Shuffle(Probability pm, double range_max)
        : Mutation(pm)
    {
        this->range_max(range_max);
    }

    void Shuffle::range_max(double rm)
    {
        if (!(0.0 <= rm && rm <= 1.0))
        {
            GA_THROW(std::invalid_argument, "The range_max parameter must be in the closed interval [0.0, 1.0].");
        }

        range_max_ = rm;
    }

    void Shuffle::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        const size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() < mutation_rate())
        {
            const size_t min_len = 2;
            const size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            const size_t range_len = rng::randomInt(min_len, max_len);

            const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
            const size_t last = first + range_len;

            std::shuffle(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last, rng::prng);
        }
    }

    Shift::Shift(Probability pm, double range_max)
        : Mutation(pm)
    {
        this->range_max(range_max);
    }

    void Shift::range_max(double rm)
    {
        if (!(0.0 <= rm && rm <= 1.0))
        {
            GA_THROW(std::invalid_argument, "The range_max parameter must be in the closed interval [0.0, 1.0].");
        }

        range_max_ = rm;
    }

    void Shift::mutate(const GA<GeneType>&, Candidate<GeneType>& candidate) const
    {
        const size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() < mutation_rate())
        {
            const size_t min_len = 2;
            const size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            const size_t range_len = rng::randomInt(min_len, max_len);

            const size_t first = rng::randomInt(0_sz, chrom_len - range_len);
            const size_t last = first + range_len;

            std::vector<GeneType> moved_elements(std::make_move_iterator(candidate.chromosome.begin()) + first,
                                                 std::make_move_iterator(candidate.chromosome.begin()) + last);
            candidate.chromosome.erase(candidate.chromosome.begin() + first,
                                       candidate.chromosome.begin() + last);

            const size_t new_pos = rng::randomInt(0_sz, candidate.chromosome.size());
            candidate.chromosome.insert(candidate.chromosome.begin() + new_pos,
                                        std::make_move_iterator(moved_elements.begin()),
                                        std::make_move_iterator(moved_elements.end()));
        }
    }

} // namespace genetic_algorithm::mutation::perm