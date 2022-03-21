/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../rng.hpp"
#include "../math.hpp"

#include <algorithm>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::real
{
    BLXa::BLXa(const std::vector<std::pair<double, double>>& bounds, double pc, double alpha) :
        BoundedCrossover(bounds, pc)
    {
        this->alpha(alpha);
    }

    void BLXa::alpha(double alpha)
    {
        if (!(0.0 <= alpha && alpha <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Alpha must be a nonnegative, finite value.");
        }

        alpha_ = alpha;
    }


    SimulatedBinary::SimulatedBinary(const std::vector<std::pair<double, double>>& bounds, double pc, double eta) :
        BoundedCrossover(bounds, pc)
    {
        this->eta(eta);
    }

    void SimulatedBinary::eta(double eta)
    {
        if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        eta_ = eta;
    }


    CandidatePair<double> Arithmetic::crossover(const GA<double>&, const Candidate<double>& parent1, const Candidate<double>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the arithmetic crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        double alpha = rng::randomReal();
        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            child1.chromosome[i] = alpha * parent1.chromosome[i] + (1.0 - alpha) * parent2.chromosome[i];
            child2.chromosome[i] = (1.0 - alpha) * parent1.chromosome[i] + alpha * parent2.chromosome[i];
        }
        /* No bounds check, the generated children's genes will always be within the bounds if the parents' genes were also within them. */

        return { child1, child2 };
    }

    CandidatePair<double> BLXa::crossover(const GA<double>&, const Candidate<double>& parent1, const Candidate<double>& parent2) const
    {
        assert(alpha_ >= 0.0);

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the BLXa crossover.");
        }
        if (parent1.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The chromosome and bounds vector lengths must be the same to perform the BLXa crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            /* Calc interval to generate the childrens genes on. */
            auto [range_min, range_max] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);
            double range_ext = alpha_ * (range_max - range_min);
            /* Generate genes from an uniform distribution on the interval. */
            child1.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
            child2.chromosome[i] = rng::randomReal(range_min - range_ext, range_max + range_ext);
            /* The children's genes might be outside the allowed interval. */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds_[i].first, bounds_[i].second);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds_[i].first, bounds_[i].second);
        }

        return { child1, child2 };
    }

    CandidatePair<double> SimulatedBinary::crossover(const GA<double>&, const Candidate<double>& parent1, const Candidate<double>& parent2) const
    {
        assert(eta_ > 0.0);

        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Simulated Binary crossover.");
        }
        if (parent1.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The chromosome and bounds vector lengths must be the same to perform the Simulated Binary crossover.");
        }

        Candidate child1{ parent1 }, child2{ parent2 };

        /* Perform crossover. */
        for (size_t i = 0; i < parent1.chromosome.size(); i++)
        {
            auto [gene_low, gene_high] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);

            /* Handle the edge case where the 2 genes are equal. */
            if ((gene_high - gene_low) <= std::numeric_limits<double>::epsilon() * std::max(std::abs(gene_low), std::abs(gene_high)))
            {
                continue;
            }

            double beta1 = 1.0 + 2.0 * (gene_low - bounds_[i].first) / (gene_high - gene_low);
            double beta2 = 1.0 + 2.0 * (bounds_[i].second - gene_high) / (gene_high - gene_low);

            double alpha1 = 2.0 - std::pow(beta1, -(eta_ + 1.0));
            double alpha2 = 2.0 - std::pow(beta2, -(eta_ + 1.0));

            double u = rng::randomReal();
            double beta1_prime = (u <= 1.0 / alpha1) ? std::pow(u * alpha1, -(eta_ + 1.0)) :
                                                       std::pow(1.0 / (2.0 - u * alpha1), -(eta_ + 1.0));

            u = rng::randomReal();
            double beta2_prime = (u <= 1.0 / alpha2) ? std::pow(u * alpha2, -(eta_ + 1.0)) :
                                                       std::pow(1.0 / (2.0 - u * alpha2), -(eta_ - 1.0));

            child1.chromosome[i] = 0.5 * (parent1.chromosome[i] + parent2.chromosome[i] - beta1_prime * std::abs(parent1.chromosome[i] - parent2.chromosome[i]));
            child2.chromosome[i] = 0.5 * (parent1.chromosome[i] + parent2.chromosome[i] + beta2_prime * std::abs(parent1.chromosome[i] - parent2.chromosome[i]));

            /* The children's genes might be outside the allowed interval. */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds_[i].first, bounds_[i].second);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds_[i].first, bounds_[i].second);
        }

        return { child1, child2 };
    }

    CandidatePair<double> Wright::crossover(const GA<double>&, const Candidate<double>& parent1, const Candidate<double>& parent2) const
    {
        if (parent1.chromosome.size() != parent2.chromosome.size())
        {
            throw std::invalid_argument("The parent chromosomes must be the same length for the Wright crossover.");
        }
        if (parent1.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The chromosome and bounds vector lengths must be the same to perform the Wright crossover.");
        }

        Candidate child1(parent1), child2(parent2);

        /* p1 is always the better parent. */
        const Candidate<double>* p1 = detail::paretoCompareLess(parent1.fitness, parent2.fitness) ? &parent2 : &parent1;
        const Candidate<double>* p2 = detail::paretoCompareLess(parent1.fitness, parent2.fitness) ? &parent1 : &parent2;
        /* Get random weights. */
        double w1 = rng::randomReal();
        double w2 = rng::randomReal();
        /* Perform crossover. */
        for (size_t i = 0; i < p1->chromosome.size(); i++)
        {
            child1.chromosome[i] = w1 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
            child2.chromosome[i] = w2 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
            /* The children's genes might be outside the allowed intervals. */
            child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds_[i].first, bounds_[i].second);
            child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds_[i].first, bounds_[i].second);
        }

        return { child1, child2 };
    }

} // namespace genetic_algorithm::crossover::real