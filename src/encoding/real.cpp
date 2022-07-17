/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../utility/rng.hpp"
#include "../selection/selection.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    RCGA::RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, const Bounds& bounds)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        this->limits(bounds);
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / this->chrom_len()));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    RCGA::RCGA(size_t chrom_len, FitnessFunction fitness_function, const Bounds& bounds)
        : RCGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function), bounds)
    {}

    RCGA::RCGA(size_t chrom_len, FitnessFunction fitness_function, const std::pair<GeneType, GeneType>& bounds)
         : RCGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function), std::vector(chrom_len, bounds))
    {}

    RCGA::RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, const std::pair<GeneType, GeneType>& bounds)
        : RCGA(pop_size, chrom_len, std::move(fitness_function), std::vector(chrom_len, bounds))
    {}

    void RCGA::limits(const Bounds& limits)
    {
        if (limits.size() != this->chrom_len())
        {
            throw std::invalid_argument("The number of limits must be equal to the chromosome length.");
        }
        if (std::any_of(limits.begin(), limits.end(), [](const std::pair<GeneType, GeneType>& b) {return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        limits_ = limits;
    }

    void RCGA::limits(const std::pair<GeneType, GeneType>& limits)
    {
        this->limits(std::vector(this->chrom_len(), limits));
    }

    const RCGA::Bounds& RCGA::limits() const& noexcept
    {
        return limits_;
    }

    RCGA::Candidate RCGA::generateCandidate() const
    {
        assert(this->chrom_len() == limits_.size());

        Candidate solution(this->chrom_len());
        for (size_t i = 0; i < solution.chromosome.size(); i++)
        {
            solution.chromosome[i] = rng::randomReal(limits_[i].first, limits_[i].second);
        }

        return solution;
    }

} // namespace genetic_algorithm