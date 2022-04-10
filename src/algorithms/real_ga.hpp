/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_REAL_GA_HPP
#define GA_REAL_GA_HPP

#include "ga_base.hpp"
#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Standard genetic algorithm that uses real encoding. \n
    * Each gene of the chromosomes is a real value.
    */
    class RCGA : public GA<double, RCGA>
    {
    public:
        /**
        * For the gene boundaries. \n
        * For example: { {gene1_min, gene1_max}, {gene2_min, gene2_max}, ... }
        */
        using Bounds = std::vector<std::pair<GeneType, GeneType>>;

        /**
        * Constructor for the RCGA.
        *
        * @param chrom_len The number of genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the real coded genes (their min and max values).
        */
        RCGA(size_t chrom_len, FitnessFunction fitnessFunction, const Bounds& bounds);

        /**
        * Constructor for the RCGA.
        *
        * @param pop_size The number of candidates in the population.
        * @param chrom_len The number of genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the real coded genes (their min and max values).
        */
        RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitnessFunction, const Bounds& bounds);

        /**
        * Sets the boundaries of the real-coded genes. \n
        * Each element must contain the lower and upper bounds of the corresponding real gene. (The min and max values of the gene.) \n
        * The number of elements must be the same as the length of the chromosomes, and the lower bounds must not be higher than the upper bounds. \n
        * Eg. in the case of chromosomes of 2 values where both genes must be between -1 and 1: \n
        * limits = { {-1.0, 1.0}, {-1.0, 1.0} }
        *
        * @param limits The boundaries of the real genes.
        */
        void limits(const Bounds& limits);

        /** @returns The current bounds set for the genes of the algorithm. */
        [[nodiscard]] const Bounds& limits() const noexcept;

    private:
        friend class GA<GeneType, RCGA>;
        Bounds limits_;

        Candidate generateCandidate() const;
    };

} // namespace genetic_algorithm

#endif // !GA_RCGA_H