/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_INTEGER_GA_HPP
#define GA_INTEGER_GA_HPP

#include "../core/ga_base.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Integer encoded genetic algorithm. \n
    * Similar to the @ref BinaryGA, but the values of the genes can be any integer in the interval [0, base), not just 0 or 1. \n
    * It also uses a slightly different mutation function with swaps and inversions.
    */
    class IntegerGA : public GA<size_t, IntegerGA>
    {
    public:
        /**
        * Construct an integer encoded genetic algorithm.
        *
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of values a gene can take. Must be > 1. If 2, same as the @ref BinaryGA.
        */
        IntegerGA(size_t chrom_len, FitnessFunction fitness_function, GeneType base);

        /**
        * Construct an integer encoded genetic algorithm.
        *
        * @param pop_size The number of candidates in a population.
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of values a gene can take. Must be > 1. If 2, same as the @ref BinaryGA.
        */
        IntegerGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, GeneType base);

        /**
        * Sets the number of values a gene can take to @p base. \n
        * The value of the base must be at least 2, and the GA is practically the same as the
        * BinaryGA if is is set to 2.
        *
        * @param base The number of values a gene can be.
        */
        void base(GeneType base);

        /** @returns The current base value set for the algorithm. */
        [[nodiscard]]
        GeneType base() const noexcept;

    private:
        friend class GA<GeneType, IntegerGA>;

        GeneType base_ = 4;

        Candidate generateCandidate() const;
    };

} // namespace genetic_algorithm

#endif // !GA_INTEGER_GA_HPP