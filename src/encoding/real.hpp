/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_REAL_HPP
#define GA_ENCODING_REAL_HPP

#include "../core/ga_base.hpp"
#include "gene_types.hpp"
#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Real coded genetic algorithm. \n
    * Each gene of the chromosomes is a floating point value.
    */
    class RCGA : public GA<RealGene>
    {
    public:
        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds are specified for each gene separately.
        *
        * @param chrom_len The number of genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the genes (their min and max values), specified for each gene.
        */
        RCGA(size_t chrom_len, FitnessFunction fitness_function, const BoundsVector& bounds);

        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds are specified for each gene separately.
        *
        * @param pop_size The number of candidates in the population.
        * @param chrom_len The number of genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the genes (their min and max values), specified for each gene.
        */
        RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, const BoundsVector& bounds);

        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds will be the same for every gene.
        *
        * @param chrom_len The number of genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries used for every gene (min, and max values).
        */
        RCGA(size_t chrom_len, FitnessFunction fitness_function, const GeneBounds& bounds);

        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds will be the same for every gene.
        *
        * @param pop_size The number of candidates in the population.
        * @param chrom_len The number of genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of every gene (min, and max values)
        */
        RCGA(size_t pop_size, size_t chrom_len, FitnessFunction fitness_function, const GeneBounds& bounds);

        /**
        * Sets the boundaries of the genes. \n
        * Each element must contain the lower and upper bounds of the corresponding gene (min and max values). \n
        * The size of @p limits must be the same as the chromosome length, and the lower bound can't be higher than the upper bound for any gene. \n
        * Eg. in the case of chromosomes of 2 values where the values of both genes should be on the closed interval [0.0, 1.0] : \n
        * limits = { {-1.0, 1.0}, {-1.0, 1.0} }
        *
        * @param limits The lower and upper boundaries of the genes.
        */
        void gene_bounds(const BoundsVector& limits);

        /**
        * Sets the the same lower and upper bounds for every gene. \n
        * The lower bound can't be higher than the upper bound. \n
        * Eg. if @p limits is { 0.0, 1.0 }, all genes will be in the closed interval [0.0, 1.0].
        *
        * @param limits The lower and upper boundaries of the genes.
        */
        void gene_bounds(const GeneBounds& limits);

    private:
        Candidate generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_REAL_HPP