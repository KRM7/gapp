/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "permutation.hpp"
#include "../core/fitness_function.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/permutation.hpp"
#include "../mutation/permutation.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <numeric>
#include <vector>
#include <memory>
#include <utility>

namespace genetic_algorithm
{
    PermutationGA::PermutationGA(std::unique_ptr<FitnessFunction<PermutationGene>> fitness_function, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        bounds_ = BoundsVector<GeneType>(this->chrom_len(), GeneBounds<GeneType>{ 0_sz, this->chrom_len() - 1 });
        crossover_method(std::make_unique<crossover::perm::Order2>());
        mutation_method(std::make_unique<mutation::perm::Inversion>(0.2));
    }

    void PermutationGA::initialize()
    {
        bounds_.resize(this->chrom_len(), GeneBounds<GeneType>{ 0_sz, this->chrom_len() - 1 });
    }

    auto PermutationGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(this->chrom_len());
        std::iota(solution.chromosome.begin(), solution.chromosome.end(), GeneType{ 0 });
        std::shuffle(solution.chromosome.begin(), solution.chromosome.end(), rng::prng);

        return solution;
    }

} // namespace genetic_algorithm