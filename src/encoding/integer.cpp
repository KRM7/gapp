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
    IntegerGA::IntegerGA(size_t chrom_len, FitnessFunction fitness_function, GeneType base)
        : IntegerGA(DEFAULT_POPSIZE, chrom_len, std::move(fitness_function), base)
    {}

    IntegerGA::IntegerGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, GeneType base)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        this->base(base);
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::integer::TwoPoint>());
        mutation_method(std::make_unique<mutation::integer::Uniform>(1.0 / this->chrom_len()));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
    }

    void IntegerGA::base(GeneType base)
    {
        if (base < 2) GA_THROW(std::invalid_argument, "The base must be at least 2.");

        base_ = base;
    }

    IntegerGA::GeneType IntegerGA::base() const noexcept
    {
        return base_;
    }

    IntegerGA::Candidate IntegerGA::generateCandidate() const
    {
        Candidate solution(this->chrom_len());
        std::generate(solution.chromosome.begin(), solution.chromosome.end(),
        [this]
        {
            return rng::randomInt<GeneType>(0, base_ - 1);
        });

        return solution;
    }

} // namespace genetic_algorithm