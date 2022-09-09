/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../utility/rng.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    RCGA::RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, const BoundsVector& bounds)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        this->gene_bounds(bounds);
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / this->chrom_len()));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    RCGA::RCGA(size_t chrom_len, FitnessFunction fitness_function, const BoundsVector& bounds)
        : RCGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function), bounds)
    {}

    RCGA::RCGA(size_t chrom_len, FitnessFunction fitness_function, const GeneBounds& bounds)
         : RCGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function), std::vector(chrom_len, bounds))
    {}

    RCGA::RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, const GeneBounds& bounds)
        : RCGA(pop_size, chrom_len, std::move(fitness_function), std::vector(chrom_len, bounds))
    {}

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

    void RCGA::initializeAlgorithmImpl()
    {
        if (bounds_.size() != chrom_len())
        {
            GA_THROW(std::logic_error, "The size of the bounds vector must match the chromosome length for the RCGA.");
        }
    }

    RCGA::Candidate RCGA::generateCandidate() const
    {
        assert(chrom_len() == bounds_.size());

        Candidate solution(chrom_len());
        for (size_t i = 0; i < solution.chromosome.size(); i++)
        {
            solution.chromosome[i] = rng::randomReal(bounds_[i].lower, bounds_[i].upper);
        }

        return solution;
    }

} // namespace genetic_algorithm