/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary.hpp"
#include "../core/fitness_function.hpp"
#include "../population/candidate.hpp"
#include "../crossover/binary.hpp"
#include "../mutation/binary.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <memory>

namespace genetic_algorithm
{
    BinaryGA::BinaryGA(std::unique_ptr<FitnessFunction<BinaryGene>> fitness_function, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        bounds_ = BoundsVector(this->chrom_len(), GeneBounds{ 0, 1 });
        crossover_method(std::make_unique<crossover::binary::TwoPoint>());
        mutation_method(std::make_unique<mutation::binary::Flip>(1.0 / this->chrom_len()));
    }

    void BinaryGA::initialize()
    {
        bounds_.resize(this->chrom_len(), GeneBounds{ 0, 1 });
    }

    auto BinaryGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(this->chrom_len());
        std::generate(solution.chromosome.begin(), solution.chromosome.end(), rng::randomBool);

        return solution;
    }

} // namespace genetic_algorithm