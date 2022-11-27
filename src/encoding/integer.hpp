/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_INTEGER_HPP
#define GA_ENCODING_INTEGER_HPP

#include "../core/ga_base.decl.hpp"
#include "../population/candidate.hpp"
#include "gene_types.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Integer encoded genetic algorithm. \n
    * Similar to the @ref BinaryGA, but the values of the genes can be any integer in the interval [offset, offset + base), not just 0 or 1. \n
    * It also uses a slightly different mutation function with swaps and inversions.
    */
    class IntegerGA : public GA<IntegerGene>
    {
    public:
        /**
        * Construct an integer encoded genetic algorithm. \n
        *
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of different values a gene can take. Must be > 1.
        * @param offset The lower bound of the genes.
        */
        IntegerGA(size_t chrom_len, FitnessFunction fitness_function, GeneType base, GeneType offset = GeneType{ 0 });

        /**
        * Construct an integer encoded genetic algorithm. \n
        *
        * @param pop_size The number of candidates in a population.
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of different values a gene can take. Must be > 1.
        * @param offset The lower bound of the genes.
        */
        IntegerGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, GeneType base, GeneType offset = GeneType{ 0 });

        /**
        * Sets the number of different values a gene can take. Must be at least 2. \n
        * The values of the genes will be integers in the range [offset, offset + base). \n
        * 
        * With base = 2 and offset = 0, the IntegerGA is effectively the same as the BinaryGA.
        *
        * @param base The number of different values allowed for a gene.
        */
        void base(GeneType base);

        /** @returns The current base set for the algorithm. */
        [[nodiscard]]
        GeneType base() const noexcept { return base_; }

        /**
        * Set an offset for the genes (the smallest integer a gene may be). \n
        * The values of the genes will be integers in the range [offset, offset + base). \n
        *
        * With base = 2 and offset = 0, the IntegerGA is effectively the same as the BinaryGA.
        * 
        * @param offset The lower bound of the genes.
        */
        void offset(GeneType offset);
        
        /** @returns The current offset set for the algorithm. */
        [[nodiscard]]
        GeneType offset() const noexcept { return offset_; }

    private:
        GeneType base_ = 4;
        GeneType offset_ = 0;

        void initialize() override;
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_INTEGER_HPP