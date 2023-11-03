/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_VARIABLE_LENGTH_HPP
#define GA_TEST_VARIABLE_LENGTH_HPP

#include "gapp.hpp"
#include <algorithm>
#include <iostream>

using namespace gapp;


class VariableCrossover final : public crossover::Crossover<BinaryGene>
{
public:
    using Crossover::Crossover;
    bool allow_variable_chrom_length() const noexcept override { return true; }
private:
    CandidatePair<GeneType> crossover(const GA<GeneType>&, const Candidate<GeneType>& parent1, const Candidate<GeneType>& parent2) const override
    {
        return { parent1, parent2 };
    }
};


class VariableMutation final : public mutation::Mutation<BinaryGene>
{
public:
    using Mutation::Mutation;
    bool allow_variable_chrom_length() const noexcept override { return true; }
private:
    void mutate(const GA<GeneType>&, const Candidate<GeneType>&, Chromosome<GeneType>& chromosome) const override
    {
        const size_t flip_count = rng::randomBinomial(chromosome.size(), mutation_rate());
        const auto flipped_indices = rng::sampleUnique(0_sz, chromosome.size(), flip_count);

        for (const auto& idx : flipped_indices) { chromosome[idx] ^= 1; }

        if (rng::randomBool()) chromosome.push_back(rng::randomBool());
    }
};


class CountOnes : public FitnessFunction<BinaryGene, 10>
{
    FitnessVector invoke(const Chromosome<BinaryGene>& x) const override
    {
        return { static_cast<double>(std::count(x.begin(), x.end(), BinaryGene{ 1 })) };
    }
};


inline void variable_chrom_length()
{
    BinaryGA ga(200);

    ga.crossover_method(VariableCrossover{ 0.8 });
    ga.mutation_method(VariableMutation{ 0.05 });

    auto solutions = ga.solve(CountOnes{}, 1000);

    std::cout << "\nVariable chromosome length problem:"
              << "\nNumber of optimal sols: " << solutions.size()
              << "\nBest fitness found: " << solutions[0].fitness[0] << "\n";
}

#endif // !GA_TEST_VARIABLE_LENGTH_HPP
