/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include "../core/fitness_function.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/integer.hpp"
#include "../mutation/integer.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <utility>
#include <memory>
#include <stdexcept>

namespace genetic_algorithm
{
    IntegerGA::IntegerGA(std::unique_ptr<FitnessFunction<IntegerGene>> fitness_function, const GeneBounds& bounds, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        this->gene_bounds(bounds);
        crossover_method(std::make_unique<crossover::integer::TwoPoint>());
        mutation_method(std::make_unique<mutation::integer::Uniform>(1.0 / this->chrom_len()));
    }

    void IntegerGA::gene_bounds(const GeneBounds& limits)
    {
        if (limits.lower > limits.upper)
        {
            GA_THROW(std::invalid_argument, "The lower bound of the genes can't be higher than the upper bound.");
        }

        bounds_ = BoundsVector(this->chrom_len(), limits);
    }

    void IntegerGA::initialize()
    {
        bounds_.resize(chrom_len(), bounds_.front());
    }

    auto IntegerGA::generateCandidate() const -> Candidate<GeneType>
    {
        Candidate<GeneType> solution(this->chrom_len());
        std::generate(solution.chromosome.begin(), solution.chromosome.end(), [this]
        {
            return rng::randomInt<GeneType>(bounds_.front().lower, bounds_.front().upper);
        });

        return solution;
    }

} // namespace genetic_algorithm