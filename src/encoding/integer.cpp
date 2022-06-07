/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../utility/rng.hpp"
#include "../selection/selection.hpp"
#include "../crossover/integer.hpp"
#include "../mutation/integer.hpp"
#include "../stop_condition/stop_condition.hpp"
#include <utility>
#include <memory>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm
{
    IntegerGA::IntegerGA(size_t chrom_len, FitnessFunction fitness_function, GeneType base)
        : IntegerGA(100, chrom_len, std::move(fitness_function), base)
    {}

    IntegerGA::IntegerGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, GeneType base)
        : GA(pop_size, chrom_len, std::move(fitness_function))
    {
        this->base(base);
        setDefaultAlgorithm();
        crossover_method(std::make_unique<crossover::integer::TwoPoint>());
        mutation_method(std::make_unique<mutation::integer::Uniform>(1.0 / chrom_len_));
        stop_condition(std::make_unique<stopping::NoEarlyStop>());
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

        Candidate solution(chrom_len_);
        std::generate(solution.chromosome.begin(), solution.chromosome.end(),
        [this]
        {
            return rng::randomInt<GeneType>(0, base_ - 1);
        });

        return solution;
    }

} // namespace genetic_algorithm