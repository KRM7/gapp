/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../encoding/real.hpp"
#include "../utility/rng.hpp"
#include "../utility/math.hpp"
#include <algorithm>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::real
{
    BLXa::BLXa(Probability pc, GeneType alpha) :
        Crossover(pc)
    {
        this->alpha(alpha);
    }

    void BLXa::alpha(GeneType alpha)
    {
        if (!(0.0 <= alpha && alpha <= std::numeric_limits<GeneType>::max()))
        {
            GA_THROW(std::invalid_argument, "Alpha must be a nonnegative, finite value.");
        }

        alpha_ = alpha;
    }


    SimulatedBinary::SimulatedBinary(Probability pc, GeneType eta) :
        Crossover(pc)
    {
        this->eta(eta);
    }

    void SimulatedBinary::eta(GeneType eta)
    {
        if (!(0.0 <= eta && eta <= std::numeric_limits<GeneType>::max()))
        {
            GA_THROW(std::invalid_argument, "Eta must be a nonnegative, finite value.");
        }

        eta_ = eta;
    }


    auto Arithmetic::crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        const auto& bounds = ga.gene_bounds();
        const size_t chrom_len = parent1.chromosome.size();

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the arithmetic crossover.");
        }
        if (bounds.size() != chrom_len)
        {
            GA_THROW(std::invalid_argument, "The chromosome and bounds vector sizes must be the same to perform the BLXa crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        const GeneType alpha = rng::randomReal();
        for (size_t i = 0; i < chrom_len; i++)
        {
            child1.chromosome[i] =    alpha      * parent1.chromosome[i] + (1.0 - alpha) * parent2.chromosome[i];
            child2.chromosome[i] = (1.0 - alpha) * parent1.chromosome[i] +     alpha     * parent2.chromosome[i];

            /* The children's genes might be outside the allowed interval (really). */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].lower, bounds[i].upper);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].lower, bounds[i].upper);
        }

        return { std::move(child1), std::move(child2) };
    }

    auto BLXa::crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        const auto& bounds = ga.gene_bounds();
        const size_t chrom_len = parent1.chromosome.size();

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the BLXa crossover.");
        }
        if (bounds.size() != chrom_len)
        {
            GA_THROW(std::invalid_argument, "The chromosome and bounds vector sizes must be the same to perform the BLXa crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < chrom_len; i++)
        {
            /* Calc interval to generate the childrens genes on. */
            const auto& [range_min, range_max] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);
            GeneType range_ext = alpha_ * (range_max - range_min);
            /* Generate genes from an uniform distribution on the interval. */
            child1.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
            child2.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
            /* The children's genes might be outside the allowed interval. */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].lower, bounds[i].upper);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].lower, bounds[i].upper);
        }

        return { std::move(child1), std::move(child2) };
    }

    auto SimulatedBinary::crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        const auto& bounds = ga.gene_bounds();
        const size_t chrom_len = parent1.chromosome.size();

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the Simulated Binary crossover.");
        }
        if (bounds.size() != chrom_len)
        {
            GA_THROW(std::invalid_argument, "The chromosome and bounds vector sizes must be the same to perform the simulated binary crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < chrom_len; i++)
        {
            const auto& [gene_low, gene_high] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);

            /* Handle the edge case where the 2 genes are equal. */
            if (math::floatIsEqual(gene_high, gene_low)) continue;

            const GeneType beta1 = 1.0 + 2.0 * (gene_low - bounds[i].lower) / (gene_high - gene_low);
            const GeneType beta2 = 1.0 + 2.0 * (bounds[i].upper - gene_high) / (gene_high - gene_low);

            const GeneType alpha1 = 2.0 - std::pow(beta1, -(eta_ + 1.0));
            const GeneType alpha2 = 2.0 - std::pow(beta2, -(eta_ + 1.0));

            const auto alphaToBetaPrime = [this](GeneType alpha)
            {
                const double u = rng::randomReal<GeneType>();

                return (u <= 1.0 / alpha) ? std::pow(u * alpha, -(eta_ + 1.0)) :
                                            std::pow(1.0 / (2.0 - u * alpha), -(eta_ + 1.0));
            };

            const GeneType beta1_prime = alphaToBetaPrime(alpha1);
            const GeneType beta2_prime = alphaToBetaPrime(alpha2);

            child1.chromosome[i] = 0.5 * (parent1.chromosome[i] + parent2.chromosome[i] - beta1_prime * std::abs(parent1.chromosome[i] - parent2.chromosome[i]));
            child2.chromosome[i] = 0.5 * (parent1.chromosome[i] + parent2.chromosome[i] + beta2_prime * std::abs(parent1.chromosome[i] - parent2.chromosome[i]));

            /* The children's genes might be outside the allowed interval. */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].lower, bounds[i].upper);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].lower, bounds[i].upper);
        }

        return { std::move(child1), std::move(child2) };
    }

    auto Wright::crossover(const GA<GeneType>& ga, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const -> CandidatePair<GeneType>
    {
        const auto& bounds = ga.gene_bounds();
        const size_t chrom_len = parent1.chromosome.size();

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            GA_THROW(std::invalid_argument, "The parent chromosomes must be the same length for the Wright crossover.");
        }
        if (bounds.size() != chrom_len)
        {
            GA_THROW(std::invalid_argument, "The chromosome and bounds vector lengths must be the same to perform the Wright crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        /* p1 is always the better parent. */
        const auto& p1 = math::paretoCompareLess(parent1.fitness, parent2.fitness) ? parent2 : parent1;
        const auto& p2 = math::paretoCompareLess(parent1.fitness, parent2.fitness) ? parent1 : parent2;

        const GeneType w1 = rng::randomReal<GeneType>();
        const GeneType w2 = rng::randomReal<GeneType>();

        for (size_t i = 0; i < chrom_len; i++)
        {
            child1.chromosome[i] = w1 * (p1.chromosome[i] - p2.chromosome[i]) + p1.chromosome[i];
            child2.chromosome[i] = w2 * (p1.chromosome[i] - p2.chromosome[i]) + p1.chromosome[i];
            /* The children's genes might be outside the allowed intervals. */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].lower, bounds[i].upper);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].lower, bounds[i].upper);
        }

        return { std::move(child1), std::move(child2) };
    }

} // namespace genetic_algorithm::crossover::real