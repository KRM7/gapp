/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "crossover_dtl.hpp"
#include "../utility/rng.hpp"
#include <stdexcept>

namespace genetic_algorithm::crossover::integer
{
    auto SinglePoint::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        return dtl::singlePointCrossoverImpl(parent1, parent2);
    }

    auto TwoPoint::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        return dtl::twoPointCrossoverImpl(parent1, parent2);
    }

    NPoint::NPoint(size_t n)
    {
        num_crossover_points(n);
    }

    void NPoint::num_crossover_points(size_t n)
    {
        if (n == 0)
        {
            throw std::invalid_argument("The number of crossover points must be at least 1 for the n-point crossover.");
        }

        n_ = n;
    }

    auto NPoint::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        if (n_ == 1)
        {
            return dtl::singlePointCrossoverImpl(parent1, parent2);
        }
        else if (n_ == 2)
        {
            return dtl::twoPointCrossoverImpl(parent1, parent2);
        }
        else
        {
            return dtl::nPointCrossoverImpl(parent1, parent2, n_);
        }
    }

    Uniform::Uniform(double ps)
    {
        swap_probability(ps);
    }

    void Uniform::swap_probability(double ps)
    {
        if (!(0.0 <= ps && ps <= 1.0))
        {
            throw std::invalid_argument("The swap probability must be in the range [0.0, 1.0] for the uniform crossover.");
        }

        ps_ = ps;
    }

    auto Uniform::crossover(const GaInfo&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        size_t chrom_len = parent1.chromosome.size();

        if (parent2.chromosome.size() != chrom_len)
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the uniform crossover.");
        }
        
        size_t num_swapped_indices = rng::randomBinomialApprox(chrom_len, ps_);
        auto swapped_indices = rng::sampleUnique(0_sz, chrom_len, num_swapped_indices);

        Candidate child1{ parent1 }, child2{ parent2 };

        for (const auto& idx : swapped_indices)
        {
            std::swap(child1.chromosome, child2.chromosome);
        }

        return { child1, child2 };
    }

} // namespace genetic_algorithm::crossover::integer