/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cstddef>

namespace genetic_algorithm::mutation::perm
{
    Inversion::Inversion(double pm, double range_max)
        : Mutation(pm)
    {
        this->range_max(range_max);
    }

    void Inversion::range_max(double rm)
    {
        if (!(0.0 <= rm && rm <= 1.0))
        {
            throw std::invalid_argument("The range_max parameter must be in the closed interval [0.0, 1.0].");
        }

        range_max_ = rm;
    }

    void Inversion::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() < pm_)
        {
            size_t min_len = 2;
            size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            size_t range_len = rng::randomInt(min_len, max_len);

            size_t first = rng::randomInt<size_t>(0, chrom_len - range_len);
            size_t last = first + range_len;

            std::reverse(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last);
        }
    }

    void Swap2::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        if (candidate.chromosome.size() < 2) return;

        if (rng::randomReal() < pm_)
        {
            auto idxs = rng::sampleUnique(candidate.chromosome.size(), 2);

            using std::swap;
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
        }
    }

    void Swap3::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        if (candidate.chromosome.size() < 3) return;

        if (rng::randomReal() < pm_)
        {
            auto idxs = rng::sampleUnique(candidate.chromosome.size(), 3);

            using std::swap;
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
            swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[2]]);
        }
    }

    Shuffle::Shuffle(double pm, double range_max)
        : Mutation(pm)
    {
        this->range_max(range_max);
    }

    void Shuffle::range_max(double rm)
    {
        if (!(0.0 <= rm && rm <= 1.0))
        {
            throw std::invalid_argument("The range_max parameter must be in the closed interval [0.0, 1.0].");
        }

        range_max_ = rm;
    }

    void Shuffle::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() < pm_)
        {
            size_t min_len = 2;
            size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            size_t range_len = rng::randomInt(min_len, max_len);

            size_t first = rng::randomInt<size_t>(0, chrom_len - range_len);
            size_t last = first + range_len;

            std::shuffle(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last, rng::prng);
        }
    }

    Shift::Shift(double pm, double range_max)
        : Mutation(pm)
    {
        this->range_max(range_max);
    }

    void Shift::range_max(double rm)
    {
        if (!(0.0 <= rm && rm <= 1.0))
        {
            throw std::invalid_argument("The range_max parameter must be in the closed interval [0.0, 1.0].");
        }

        range_max_ = rm;
    }

    void Shift::mutate(const GaInfo&, Candidate<GeneType>& candidate) const
    {
        size_t chrom_len = candidate.chromosome.size();

        if (chrom_len < 2) return;

        if (rng::randomReal() < pm_)
        {
            size_t min_len = 2;
            size_t max_len = std::max(size_t(range_max_ * chrom_len), min_len);
            size_t range_len = rng::randomInt(min_len, max_len);
            size_t first = rng::randomInt<size_t>(0, chrom_len - range_len);
            size_t last = first + range_len;

            std::vector<GeneType> moved_elements(std::make_move_iterator(candidate.chromosome.begin()) + first,
                                                 std::make_move_iterator(candidate.chromosome.begin()) + last);
            candidate.chromosome.erase(candidate.chromosome.begin() + first,
                                       candidate.chromosome.begin() + last);

            size_t new_pos = rng::randomInt(size_t{ 0 }, candidate.chromosome.size());
            candidate.chromosome.insert(candidate.chromosome.begin() + new_pos,
                                        std::make_move_iterator(moved_elements.begin()),
                                        std::make_move_iterator(moved_elements.end()));
        }
    }

} // namespace genetic_algorithm::mutation::perm