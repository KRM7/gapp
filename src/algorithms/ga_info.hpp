/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_GA_INFO_HPP
#define GA_GA_INFO_HPP

#include <vector>
#include <atomic>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo
    {
    public:
        GaInfo(size_t chrom_len);

        GaInfo(const GaInfo&) = default;
        GaInfo(GaInfo&&) = default;
        GaInfo& operator=(const GaInfo&) = default;
        GaInfo& operator=(GaInfo&&) = default;
        virtual ~GaInfo() = default;

        /**
        * Should be set to false if the fitness function does not change over time. \n
        * (The fitness function will always return the same fitness values for a given chromosome.) \n
        * Used to eliminate unnecesary objective function evaluations.
        */
        bool dynamic_fitness = false;

        /**
        * All pareto optimal optimal solutions found by the algorithm will be stored in the solutions,
        * not just the ones in the current population if set to true. \n
        */
        bool archive_optimal_solutions = false;

        /**
        * Sets the length of the chromosomes (number of genes) of the Candidate solutions used in the algorithm to @p len. \n
        * The chromosome length must be at least 1.
        *
        * @param len The length of the chromosomes.
        */
        void chrom_len(size_t len);

        /** @returns The chromosome length used for the candidates of the population. */
        [[nodiscard]] size_t chrom_len() const noexcept;

        /**
        * Sets the number of Candidates used in the population to @p size. \n
        * The population size must be at least 1.
        *
        * @param size The size of the populations.
        */
        void population_size(size_t size);

        /** @returns The population size of the algorithm. */
        [[nodiscard]] size_t population_size() const noexcept;

        /** @returns The maximum number of generations set for the algorithm. */
        [[nodiscard]] size_t max_gen() const noexcept;

        /** @returns The number of objectives. Determined by the algorithm and returns 0 before the start of the algorithm. */
        [[nodiscard]] size_t num_objectives() const noexcept;

        /** @returns The fitness matrix of the population. */
        [[nodiscard]]
        virtual std::vector<std::vector<double>> fitness_matrix() const = 0;

        /** @returns The number of fitness evaluations performed so far by the algorithm. */
        [[nodiscard]] size_t num_fitness_evals() const noexcept;

        /** @returns The current value of the generation counter. */
        [[nodiscard]] size_t generation_cntr() const noexcept;

    protected:

        size_t generation_cntr_ = 0;
        size_t num_objectives_ = 0;
        std::atomic<size_t> num_fitness_evals_ = 0;

        size_t chrom_len_;
        size_t population_size_ = 100;
        size_t max_gen_ = 500;

        void max_gen(size_t max_gen);
        void num_objectives(size_t n);
    };

} // namespace genetic_algorithm

#endif // !GA_GA_INFO_HPP