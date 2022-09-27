/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_BINARY_HPP
#define GA_ENCODING_BINARY_HPP

#include "../core/ga_base.hpp"
#include "gene_types.hpp"
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
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function to find the maximum of in the algorithm.
        */
        BinaryGA(size_t chrom_len, FitnessFunction fitness_function);

        /**
        * Construct a binary encoded genetic algorithm.
        *
        * @param pop_size The number of candidates in the population.
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function to find the maximum of in the algorithm.
        */
        BinaryGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function);

    private:
        void initialize() override;
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_BINARY_HPP