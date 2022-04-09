/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cstddef>

namespace genetic_algorithm::mutation::perm
{
    void Inversion::mutate(const GaInfo&, Candidate<size_t>& candidate) const
    {
        if (candidate.chromosome.size() > 1 && rng::randomReal() < pm_)
        {
            size_t len = rng::randomInt(size_t{ 2 }, std::max(size_t(0.75 * candidate.chromosome.size()), size_t{ 2 }));
            size_t first = rng::randomInt<size_t>(0, candidate.chromosome.size() - len);
            size_t last = first + len;

            std::reverse(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last);
        }
    }

    void Swap2::mutate(const GaInfo&, Candidate<size_t>& candidate) const
    {
        if (candidate.chromosome.size() < 2)
        {
            throw std::invalid_argument("The chromosome must have at least 2 genes for the Swap2 mutation.");
        }

        if (rng::randomReal() < pm_)
        {
            auto idxs = rng::sampleUnique(candidate.chromosome.size(), 2);
            std::swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
        }
    }

    void Swap3::mutate(const GaInfo&, Candidate<size_t>& candidate) const
    {
        if (candidate.chromosome.size() < 3)
        {
            throw std::invalid_argument("The chromosome must have at least 3 genes for the Swap3 mutation.");
        }

        if (rng::randomReal() < pm_)
        {
            auto idxs = rng::sampleUnique(candidate.chromosome.size(), 3);
            std::swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[1]]);
            std::swap(candidate.chromosome[idxs[0]], candidate.chromosome[idxs[2]]);
        }
    }

    void Shuffle::mutate(const GaInfo&, Candidate<size_t>& candidate) const
    {
        if (candidate.chromosome.size() > 1 && rng::randomReal() < pm_)
        {
            size_t len = rng::randomInt(size_t{ 2 }, std::max(size_t(0.5 * candidate.chromosome.size()), size_t{ 2 }));
            size_t first = rng::randomInt<size_t>(0, candidate.chromosome.size() - len);
            size_t last = first + len;

            std::shuffle(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last, rng::prng);
        }
    }

    void Shift::mutate(const GaInfo&, Candidate<size_t>& candidate) const
    {
        if (candidate.chromosome.size() > 1 && rng::randomReal() < pm_)
        {
            size_t len = rng::randomInt(size_t{ 2 }, std::max(size_t(0.75 * candidate.chromosome.size()), size_t{ 2 }));
            size_t first = rng::randomInt<size_t>(0, candidate.chromosome.size() - len);
            size_t last = first + len;

            std::vector<size_t> moved_elements(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last);
            candidate.chromosome.erase(candidate.chromosome.begin() + first, candidate.chromosome.begin() + last);

            size_t new_pos = rng::randomInt(size_t{ 0 }, candidate.chromosome.size());
            candidate.chromosome.insert(candidate.chromosome.begin() + new_pos, moved_elements.begin(), moved_elements.end());
        }
    }

} // namespace genetic_algorithm::mutation::perm