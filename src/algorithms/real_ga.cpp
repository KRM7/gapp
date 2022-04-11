/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real_ga.hpp"
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
    RCGA::RCGA(size_t chrom_len, FitnessFunction fitnessFunction, const Bounds& bounds)
        : GA(chrom_len, std::move(fitnessFunction))
    {
        this->limits(bounds);

        num_objectives(getNumObjectives(fitness_function_));

        if (num_objectives() == 1)
        {
            selection_method(std::make_unique<selection::single_objective::Tournament>());
        }
        else
        {
            selection_method(std::make_unique<selection::multi_objective::NSGA3>());
        }
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / chrom_len));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    RCGA::RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitnessFunction, const Bounds& bounds)
        : GA(pop_size, chrom_len, std::move(fitnessFunction))
    {
        this->limits(bounds);

        num_objectives(getNumObjectives(fitness_function_));

        if (num_objectives() == 1)
        {
            selection_method(std::make_unique<selection::single_objective::Tournament>());
        }
        else
        {
            selection_method(std::make_unique<selection::multi_objective::NSGA3>());
        }
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / chrom_len));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    void RCGA::limits(const Bounds& limits)
    {
        if (limits.size() != chrom_len_)
        {
            throw std::invalid_argument("The number of limits must be equal to the chromosome length.");
        }
        if (std::any_of(limits.begin(), limits.end(), [](std::pair<GeneType, GeneType> b) {return b.first > b.second; }))
        {
            throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
        }

        limits_ = limits;
    }

    const RCGA::Bounds& RCGA::limits() const noexcept
    {
        return limits_;
    }

    RCGA::Candidate RCGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);
        assert(chrom_len_ == limits_.size());

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomReal(limits_[i].first, limits_[i].second));
        }

        return sol;
    }

} // namespace genetic_algorithm