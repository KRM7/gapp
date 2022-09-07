/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../utility/rng.hpp"
#include "../crossover/permutation.hpp"
#include "../mutation/permutation.hpp"
#include "../stop_condition/stop_condition.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <memory>
#include <stdexcept>

namespace genetic_algorithm
{
    PermutationGA::PermutationGA(size_t chrom_len, FitnessFunction fitness_function)
        : PermutationGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function))
    {}

    PermutationGA::PermutationGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        bounds_ = BoundsVector(chrom_len, GeneBounds(0, chrom_len - 1));
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::perm::Order2>());
        mutation_method(std::make_unique<mutation::perm::Inversion>(0.2));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    PermutationGA::Candidate PermutationGA::generateCandidate() const
    {
        Candidate solution(this->chrom_len());
        std::iota(solution.chromosome.begin(), solution.chromosome.end(), GeneType{ 0 });
        std::shuffle(solution.chromosome.begin(), solution.chromosome.end(), rng::prng);

        return solution;
    }

} // namespace genetic_algorithm