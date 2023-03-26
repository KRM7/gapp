/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_REAL_HPP
#define GA_ENCODING_REAL_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../core/fitness_function.hpp"
#include "../population/candidate.hpp"
#include "../crossover/real.hpp"
#include "../mutation/real.hpp"
#include "../utility/bounded_value.hpp"
#include <concepts>
#include <memory>
#include <utility>
#include <cstddef>

namespace genetic_algorithm
{
    /**
    * Real coded genetic algorithm. \n
    * Each gene of the chromosomes is a floating point value.
    */
    class RCGA final : public GA<RealGene>
    {
    public:
        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds are specified for each gene separately.
        *
        * @param fitness_function The fitness function used in the algorithm. Can't be a nullptr.
        * @param bounds The boundaries of the genes (their min and max values), specified for each gene.
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        RCGA(std::unique_ptr<FitnessFunction<RealGene>> fitness_function, BoundsVector<GeneType> bounds, Positive<size_t> population_size = DEFAULT_POPSIZE);

        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds will be the same for every gene.
        *
        * @param fitness_function The fitness function used in the algorithm. Can't be a nullptr.
        * @param bounds The boundaries of the genes (their min and max values).
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        RCGA(std::unique_ptr<FitnessFunction<RealGene>> fitness_function, GeneBounds<GeneType> bounds, Positive<size_t> population_size = DEFAULT_POPSIZE);

        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds are specified for each gene separately.
        *
        * @param fitness_function The fitness function used in the algorithm.
        * @param bounds The boundaries of the genes (their min and max values), specified for each gene.
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        template<typename F>
        requires std::derived_from<F, FitnessFunction<GeneType>> && std::is_final_v<F>
        RCGA(F fitness_function, BoundsVector<GeneType> bounds, Positive<size_t> population_size = DEFAULT_POPSIZE) :
            GA(std::make_unique<F>(std::move(fitness_function)), population_size)
        {
            gene_bounds(std::move(bounds));
            crossover_method(std::make_unique<crossover::real::Wright>());
            mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / chrom_len()));
        }

        /**
        * Construct a real encoded genetic algorithm. \n
        * The gene bounds will be the same for every gene.
        *
        * @param fitness_function The fitness function used in the algorithm.
        * @param bounds The boundaries of the genes (their min and max values).
        * @param population_size The number of candidates in the population. Must be at least 1.
        */
        template<typename F>
        requires std::derived_from<F, FitnessFunction<GeneType>> && std::is_final_v<F>
        RCGA(F fitness_function, GeneBounds<GeneType> bounds, Positive<size_t> population_size = DEFAULT_POPSIZE) :
            GA(std::make_unique<F>(std::move(fitness_function)), population_size)
        {
            gene_bounds(bounds);
            crossover_method(std::make_unique<crossover::real::Wright>());
            mutation_method(std::make_unique<mutation::real::Gauss>(1.0 / chrom_len()));
        }


        /**
        * Sets the boundaries of the genes. \n
        * Each element must contain the lower and upper bounds of the corresponding gene (min and max values). \n
        * The size of @p bounds must be the same as the chromosome length, and the lower bound can't be higher than the upper bound for any gene. \n
        * Eg. in the case of chromosomes of 2 values where the values of both genes should be on the closed interval [0.0, 1.0] : \n
        * bounds = { {-1.0, 1.0}, {-1.0, 1.0} }
        *
        * @param bounds The lower and upper boundaries of the genes.
        */
        void gene_bounds(BoundsVector<GeneType> bounds);

        /**
        * Sets the the same lower and upper bounds for every gene. \n
        * The lower bound can't be higher than the upper bound. \n
        * Eg. if @p bounds is { 0.0, 1.0 }, all genes will be in the closed interval [0.0, 1.0].
        *
        * @param bounds The lower and upper boundaries of the genes.
        */
        void gene_bounds(GeneBounds<GeneType> bounds);

        using GA::gene_bounds;

    private:
        void initialize() override;
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_REAL_HPP