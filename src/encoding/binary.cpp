/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../utility/rng.hpp"
#include "../crossover/binary.hpp"
#include "../mutation/binary.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <utility>
#include <memory>
#include <cassert>

namespace genetic_algorithm
{
    BinaryGA::BinaryGA(size_t chrom_len, FitnessFunction fitness_function)
        : BinaryGA(100, chrom_len, std::move(fitness_function))
    {}

    BinaryGA::BinaryGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::binary::TwoPoint>());
        mutation_method(std::make_unique<mutation::binary::Flip>(1.0 / chrom_len_));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    BinaryGA::Candidate BinaryGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);

        Candidate solution(chrom_len_);
        std::generate(solution.chromosome.begin(), solution.chromosome.end(), rng::randomBool);

        return solution;
    }

} // namespace genetic_algorithm