/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_PERMUTATION_HPP
#define GA_ENCODING_PERMUTATION_HPP

#include "../core/ga_base.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Genetic algorithm that uses permutational encoding. \n
    * The genes of the chromosomes are all unique unsigned integers on thhe closed interval [0, chrom_len - 1].
    */
    class PermutationGA : public GA<size_t, PermutationGA>
    {
    public:
        /**
        * Construct a permutation encoded genetic algorithm.
        *
        * @param chrom_len The number of genes in the chromosomes.
        * @param fitness_function The fitness function used in the algorithm to find the maximum of.
        */
        PermutationGA(size_t chrom_len, FitnessFunction fitness_function);

        /**
        * Construct a permutation encoded genetic algorithm.
        *
        * @param pop_size The number of candidates in the population.
        * @param chrom_len The number of genes in the chromosomes.
        * @param fitness_function The fitness function used in the algorithm to find the maximum of.
        */
        PermutationGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function);

    private:
        friend class GA<GeneType, PermutationGA>;

        Candidate generateCandidate() const;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_PERMUTATION_HPP