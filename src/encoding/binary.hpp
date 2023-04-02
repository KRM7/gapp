/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_BINARY_HPP
#define GA_ENCODING_BINARY_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/binary.hpp"
#include "../mutation/binary.hpp"
#include "../utility/bounded_value.hpp"
#include <concepts>
#include <memory>
#include <utility>
#include <cstddef>

namespace genetic_algorithm
{
    /** Standard binary encoded genetic algorithm. */
    class BinaryGA final : public GA<BinaryGene>
    {
    public:
        /**
        * Construct a binary encoded genetic algorithm.
        *
        * @param fitness_function The fitness function used in the algorithm. Can't be a nullptr.
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        explicit BinaryGA(std::unique_ptr<FitnessFunctionBase<BinaryGene>> fitness_function, Positive<size_t> population_size = DEFAULT_POPSIZE);

        /**
        * Construct a binary encoded genetic algorithm.
        *
        * @param fitness_function The fitness function used in the algorithm.
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        template<typename F>
        requires std::derived_from<F, FitnessFunctionBase<GeneType>> && std::is_final_v<F>
        explicit BinaryGA(F fitness_function, Positive<size_t> population_size = DEFAULT_POPSIZE) :
            GA(std::make_unique<F>(std::move(fitness_function)), population_size)
        {
            bounds_ = BoundsVector<GeneType>(chrom_len(), GeneBounds<GeneType>{ 0, 1 });
            crossover_method(std::make_unique<crossover::binary::TwoPoint>());
            mutation_method(std::make_unique<mutation::binary::Flip>(1.0 / chrom_len()));
        }

    private:
        void initialize() override;
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_BINARY_HPP