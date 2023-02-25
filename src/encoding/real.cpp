/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../core/fitness_function.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <memory>
#include <stdexcept>
#include <utility>
#include <cassert>

namespace genetic_algorithm
{
    RCGA::RCGA(std::unique_ptr<FitnessFunction<RealGene>> fitness_function, const BoundsVector& bounds, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        this->gene_bounds(bounds);
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / this->chrom_len()));
    }

    RCGA::RCGA(std::unique_ptr<FitnessFunction<RealGene>> fitness_function, const GeneBounds& bounds, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        this->gene_bounds(bounds);
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / this->chrom_len()));
    }


    void RCGA::gene_bounds(const BoundsVector& limits)
    {
        if (limits.size() != this->chrom_len())
        {
            GA_THROW(std::invalid_argument, "The number of limits must be equal to the chromosome length.");
        }
        if (std::any_of(limits.begin(), limits.end(), [](const auto& limit) { return limit.lower > limit.upper; }))
        {
            GA_THROW(std::invalid_argument, "The lower bound must be lower than the upper bound for each gene.");
        }

        bounds_ = limits;
    }

    void RCGA::gene_bounds(const GeneBounds& limits)
    {
        gene_bounds(BoundsVector(chrom_len(), limits));
    }

    void RCGA::initialize()
    {
        if (bounds_.size() != chrom_len())
        {
            GA_THROW(std::logic_error, "The size of the bounds vector must match the chromosome length for the RCGA.");
        }
    }

    auto RCGA::generateCandidate() const -> Candidate<GeneType>
    {
        assert(chrom_len() == bounds_.size());

        Candidate<GeneType> solution(chrom_len());
        for (size_t i = 0; i < solution.chromosome.size(); i++)
        {
            solution.chromosome[i] = rng::randomReal(bounds_[i].lower, bounds_[i].upper);
        }

        return solution;
    }

} // namespace genetic_algorithm