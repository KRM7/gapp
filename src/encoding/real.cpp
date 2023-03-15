/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../core/fitness_function.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <vector>
#include <memory>
#include <utility>

namespace genetic_algorithm
{
    RCGA::RCGA(std::unique_ptr<FitnessFunction<RealGene>> fitness_function, const BoundsVector<GeneType>& bounds, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        this->gene_bounds(bounds);
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / this->chrom_len()));
    }

    RCGA::RCGA(std::unique_ptr<FitnessFunction<RealGene>> fitness_function, GeneBounds<GeneType> bounds, size_t population_size) :
        GA(std::move(fitness_function), population_size)
    {
        this->gene_bounds(bounds);
        crossover_method(std::make_unique<crossover::real::Wright>());
        mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / this->chrom_len()));
    }

    void RCGA::gene_bounds(const BoundsVector<GeneType>& bounds)
    {
        GA_ASSERT(bounds.size() == this->chrom_len(), "The size of the bounds vector must match the chromosome length.");

        bounds_ = bounds;
    }

    void RCGA::gene_bounds(GeneBounds<GeneType> bounds)
    {
        bounds_ = BoundsVector<GeneType>(chrom_len(), bounds);
    }

    void RCGA::initialize()
    {
        GA_ASSERT(bounds_.size() == chrom_len(), "The size of the bounds vector must match the chromosome length.");
    }

    auto RCGA::generateCandidate() const -> Candidate<GeneType>
    {
        GA_ASSERT(chrom_len() == bounds_.size(), "The size of the bounds vector must match the chromosome length.");

        Candidate<GeneType> solution(chrom_len());
        for (size_t i = 0; i < solution.chromosome.size(); i++)
        {
            solution.chromosome[i] = rng::randomReal(bounds_[i].lower(), bounds_[i].upper());
        }

        return solution;
    }

} // namespace genetic_algorithm