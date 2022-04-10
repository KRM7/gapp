/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer_ga.hpp"
#include "../utility/rng.hpp"
#include "../selection/selection.hpp"
#include "../crossover/integer.hpp"
#include "../mutation/integer.hpp"
#include <utility>
#include <memory>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    IntegerGA::IntegerGA(size_t chrom_len, FitnessFunction fitnessFunction, GeneType base)
        : GA(chrom_len, std::move(fitnessFunction))
    {
        this->base(base);

        num_objectives(getNumObjectives(fitness_function_));

        if (num_objectives() == 1)
        {
            selection_method(std::make_unique<selection::single_objective::Tournament>());
        }
        else
        {
            selection_method(std::make_unique<selection::multi_objective::NSGA3>());
        }
        crossover_method(std::make_unique<crossover::integer::TwoPoint>());
        mutation_method(std::make_unique<mutation::integer::Uniform>(1.0 / chrom_len));
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