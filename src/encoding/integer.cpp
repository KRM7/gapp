/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include "../crossover/integer.hpp"
#include "../mutation/integer.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <utility>
#include <memory>
#include <stdexcept>

namespace genetic_algorithm
{
    IntegerGA::IntegerGA(size_t chrom_len, FitnessFunction fitness_function, GeneType base, GeneType offset)
        : IntegerGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function), base, offset)
    {}

    IntegerGA::IntegerGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, GeneType base, GeneType offset)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        bounds_ = BoundsVector(chrom_len, GeneBounds(offset, offset + base - 1));
        this->base(base);
        this->offset(offset);
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::integer::TwoPoint>());
        mutation_method(std::make_unique<mutation::integer::Uniform>(1.0 / this->chrom_len()));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    void IntegerGA::base(GeneType base)
    {
        if (base < 2) GA_THROW(std::invalid_argument, "The base must be at least 2.");

        base_ = base;
        bounds_ = BoundsVector(chrom_len(), GeneBounds(offset_, offset_ + base - 1));
    }

    void IntegerGA::offset(GeneType offset)
    {
        offset_ = offset;
        bounds_ = BoundsVector(chrom_len(), GeneBounds(offset, offset + base_ - 1));
    }

    void IntegerGA::initialize()
    {
        bounds_.resize(chrom_len(), GeneBounds(offset_, offset_ + base_ - 1));
    }

    auto IntegerGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(this->chrom_len());
        std::generate(solution.chromosome.begin(), solution.chromosome.end(), [this]
        {
            return rng::randomInt<GeneType>(offset_, offset_ + base_ - 1);
        });

        return solution;
    }

} // namespace genetic_algorithm