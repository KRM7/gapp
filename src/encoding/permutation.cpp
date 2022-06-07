/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../utility/rng.hpp"
#include "../selection/selection.hpp"
#include "../crossover/permutation.hpp"
#include "../mutation/permutation.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <cassert>

namespace genetic_algorithm
{
    void PermutationGA::setDefaultOperators()
    {
        if (num_objectives() == 1)
        {
            selection_method(std::make_unique<selection::single_objective::Tournament>());
        }
        else
        {
            selection_method(std::make_unique<selection::multi_objective::NSGA3>());
        }
        crossover_method(std::make_unique<crossover::perm::Order2>());
        mutation_method(std::make_unique<mutation::perm::Inversion>(0.2));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    PermutationGA::PermutationGA(size_t chrom_len, FitnessFunction fitnessFunction)
        : GA(chrom_len, std::move(fitnessFunction))
    {
        num_objectives(getNumObjectives(fitness_function_));
        setDefaultOperators();
    }

    PermutationGA::PermutationGA(size_t pop_size, size_t chrom_len, FitnessFunction fitnessFunction)
        : GA(pop_size, chrom_len, std::move(fitnessFunction))
    {
        num_objectives(getNumObjectives(fitness_function_));
        setDefaultOperators();
    }

    PermutationGA::Candidate PermutationGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate sol;

        Chromosome chrom(chrom_len_);
        std::iota(chrom.begin(), chrom.end(), GeneType{ 0 });
        std::shuffle(chrom.begin(), chrom.end(), rng::prng);

        sol.chromosome = chrom;

        return sol;
    }

} // namespace genetic_algorithm