/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_BINARY_HPP
#define GA_ENCODING_BINARY_HPP

#include "../core/ga_base.decl.hpp"
#include "../population/candidate.hpp"
#include "gene_types.hpp"
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
        * @param fitness_function The fitness function used in the algorithm.
        * @param population_size The number of candidates in the population.
        */
        explicit BinaryGA(std::unique_ptr<FitnessFunction<BinaryGene>> fitness_function, size_t population_size = DEFAULT_POPSIZE);

        /**
        * Construct a binary encoded genetic algorithm.
        *
        * @param fitness_function The fitness function used in the algorithm.
        * @param population_size The number of candidates in the population.
        */
        template<typename F>
        requires FitnessFunctionType<F, BinaryGene> && std::is_final_v<F>
        explicit BinaryGA(F fitness_function, size_t population_size = DEFAULT_POPSIZE) :
            GA(std::make_unique<F>(std::move(fitness_function)), population_size)
        {}

    private:
        void initialize() override;
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_BINARY_HPP