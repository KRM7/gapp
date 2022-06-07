/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BINARY_GA_HPP
#define GA_BINARY_GA_HPP

#include "../core/ga_base.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Standard genetic algorithm with binary encoding. \n
    * The genes are encoded as char.
    */
    class BinaryGA final : public GA<char, BinaryGA>
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
        friend class GA<GeneType, BinaryGA>;

        Candidate generateCandidate() const;
    };

} // namespace genetic_algorithm

#endif // !GA_BINARY_GA_HPP