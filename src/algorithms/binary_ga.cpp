/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary_ga.hpp"
#include "../utility/rng.hpp"
#include "../selection/selection.hpp"
#include "../crossover/binary.hpp"
#include "../mutation/binary.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <utility>
#include <memory>
#include <cassert>
#include <cstdlib>

namespace genetic_algorithm
{
    void BinaryGA::setDefaultOperators()
    {
        if (num_objectives() == 1)
        {
            selection_method(std::make_unique<selection::single_objective::Tournament>());
        }
        else
        {
            selection_method(std::make_unique<selection::multi_objective::NSGA3>());
        }
        crossover_method(std::make_unique<crossover::binary::TwoPoint>());
        mutation_method(std::make_unique<mutation::binary::Flip>(1.0 / chrom_len_));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    BinaryGA::BinaryGA(size_t chrom_len, FitnessFunction fitness_function)
        : GA(chrom_len, std::move(fitness_function))
    {
        num_objectives(getNumObjectives(fitness_function_));
        setDefaultOperators();
    }

    BinaryGA::BinaryGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        num_objectives(getNumObjectives(fitness_function_));
        setDefaultOperators();
    }

    BinaryGA::Candidate BinaryGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(GeneType{ rng::randomBool() });
        }

        return sol;
    }

} // namespace genetic_algorithm