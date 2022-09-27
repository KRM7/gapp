/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_PERMUTATION_HPP
#define GA_ENCODING_PERMUTATION_HPP

#include "../core/ga_base.hpp"
#include "gene_types.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Genetic algorithm where the chromosomes encode permutations. \n
    * The genes of the chromosomes are all unique unsigned integers on thhe closed interval [0, chrom_len - 1].
    */
    class PermutationGA : public GA<PermutationGene>
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
        void initialize() override;
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_PERMUTATION_HPP