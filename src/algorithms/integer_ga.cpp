/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer_ga.hpp"
#include "../utility/rng.hpp"
#include <utility>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    IntegerGA::IntegerGA(size_t chrom_len, FitnessFunction fitnessFunction, GeneType base)
        : GA(chrom_len, std::move(fitnessFunction))
    {
        this->base(base);
    }

    void IntegerGA::base(GeneType base)
    {
        if (base < 2)
        {
            throw std::invalid_argument("The base must be at least 2.");
        }
        base_ = base;
    }

    IntegerGA::GeneType IntegerGA::base() const noexcept
    {
        return base_;
    }

    IntegerGA::Candidate IntegerGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomInt(GeneType{ 0 }, base_ - 1));
        }

        return sol;
    }

} // namespace genetic_algorithm