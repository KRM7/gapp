/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this softwareand associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright noticeand this permission notice shall be included in all
*  copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
*/

/**
* This file contains the base genetic algorithm class template.
*
* @file base_ga.h
*/

#ifndef GA_BASE_GA_H
#define GA_BASE_GA_H

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <utility>
#include <functional>
#include <atomic>
#include <cstddef>

/** Genetic algorithms and random number generation. */
namespace genetic_algorithm
{
    /**
    * Base GA class.
    *
    * @tparam geneType The type of the genes in the candidates' chromosomes.
    */
    template<typename geneType>
    class GA
    {
    public:

        /** Structure containing stats of the single-objective algorithm. */
        struct History
        {
            std::vector<double> fitness_mean;	/**< The mean fitness values of each generation. */
            std::vector<double> fitness_sd;		/**< The standard deviation of the fitness values of each generation. */
            std::vector<double> fitness_min;	/**< The lowest fitness value in each generation. */
            std::vector<double> fitness_max;	/**< The highest fitness value in each generation. */

            void clear() noexcept;
            void reserve(size_t new_capacity);
            void add(double mean, double sd, double min, double max);
        };

        /** The candidates used in the algorithm, each representing a solution to the problem. */
        struct Candidate
        {
            std::vector<geneType> chromosome;	/**< The chromosome encoding the solution. */
            std::vector<double> fitness;		/**< The fitness values of the candidate solution. */

            double selection_pdf = 0.0;			/**< The probability of selecting the candidate (SOGA). */
            double selection_cdf = 0.0;			/**< The value of the cumulative distribution function for the candidate (SOGA). */

            size_t rank = 0;					/**< Non-domination rank (used in both the NSGA-II and NSGA-III). */
            double distance = 0.0;				/**< Crowding distance (NSGA-II), or the distance to closest reference point (NSGA-III). */
            size_t ref_idx = 0;					/**< Index of the associated reference point (NSGA-III). */
            size_t niche_count = 0;				/**< Number of candidates associated with the same reference point as this candidate (NSGA-III). */

            bool is_evaluated = false;			/**< False if the candidate's fitness value needs to be computed. */

            Candidate();
            Candidate(const std::vector<geneType>& chrom);
            Candidate(std::vector<geneType>&& chrom) noexcept;

            bool operator==(const Candidate& rhs) const;
            bool operator!=(const Candidate& rhs) const;
        };

        /**
        * Hasher for the candidates so they can be stored in an unordered set. \n
        * A hash function must be defined for the geneType for this to work.
        */
        struct CandidateHasher
        {
            size_t operator()(const Candidate& c) const noexcept;
        };

        using Chromosome = std::vector<geneType>;								/**< . */
        using CandidatePair = std::pair<Candidate, Candidate>;					/**< . */
        using CandidateVec = std::vector<Candidate>;							/**< . */
        using CandidateSet = std::unordered_set<Candidate, CandidateHasher>;	/**< . */
        using Population = std::vector<Candidate>;								/**< . */

        using fitnessFunction_t = std::function<std::vector<double>(const Chromosome&)>;	/**< The type of the fitness function. */
        using selectionFunction_t = std::function<Candidate(const Population&)>;			/**< The type of the selection function. */
        using crossoverFunction_t = std::function<CandidatePair(const Candidate&, const Candidate&, double)>;	/**< The type of the crossover function. */
        using mutationFunction_t = std::function<void(Candidate&, double)>;					/**< The type of the mutation function. */
        using repairFunction_t = std::function<Chromosome(const Chromosome&)>;				/**< The type of the repair function. */
        using callbackFunction_t = std::function<void(const GA*)>;

        /**
        * The type of the genetic algorithm used, depending on the problem type (single- or multi-objective optimization). \n
        * Set the mode used with @ref mode.
        */
        enum class Mode
        {
            single_objective,			/**< Simple single-objective genetic algorithm. */
            multi_objective_sorting,	/**< Non-dominated sorting genetic algorithm (NSGA-II) for multi-objective optimization. */
            multi_objective_decomp		/**< NSGA-III algorithm for many-objective optimization. */
        };

        /**
        * The possible stop conditions used in the algorithm. The algorithm always stops when @ref max_gen has been reached,
        * regardless of the stop condition selected. \n
        * Some of the stop condition do not work for multi-objective problems (fitness_mean_stall and fitness_best_stall).
        * Choose the stop condition with @ref stop_condition. \n
        */
        enum class StopCondition
        {
            max_gen,			/**< Only stop when @ref max_gen is reached. */
            fitness_value,		/**< Stop when a solution was found which dominates a reference fitness value. @see fitness_threshold */
            fitness_evals,		/**< Stop when the fitness function has been evaluated a set number of times. @see max_fitness_evals */
            fitness_mean_stall,	/**< Stop when the mean fitness of the population doesn't improve at least @ref stall_threshold over @ref stall_gen_count. */
            fitness_best_stall	/**< Stop when the highest fitness of the population doesn't improve at least @ref stall_threshold over @ref stall_gen_count. */
        };

        /**
        * The possible selection methods used in the single-objective algorithm. \n
        * Choose the selection method with @ref selection_method. \n
        * If the custom option is selected, the function set with @ref setWeightCalcFunction will be used to calculate
        * the selection weights of the candidates, and weight proportional selection will be used to select candidates. \n
        * All of the methods work with negative fitness values.
        */
        enum class SogaSelection
        {
            roulette,		/**< Standard roulette selection adapted to also work with negative fitness values. No parameters. */
            rank,			/**< Standard rank selection. @see rank_sel_min_w_, @see rank_sel_max_w_ */
            tournament,		/**< Standard tournament selection. @see tournament_size */
            sigma,			/**< Sigma fitness scaling. @see sigma_scale */
            boltzmann,		/**< Standard Boltzmann selection. @see boltzmann_temps */
            custom			/**< A user defined function is used to compute the selection probabilities. @see customCalcWeights */
        };

        /**
        * Should be set to false if the fitness function does not change over time. \n
        * (The fitness function will always return the same value for a given chromosome.) \n
        * Used to eliminate unnecesary fitness evaluations.
        */
        bool changing_fitness_func = false;

        /**
        * All pareto optimal optimal solutions found in the algorithm will be stored in the solutions,
        * not just the ones in the current population if this is set to true. \n
        * Setting it to false can speed up the algorithm.
        */
        bool archive_optimal_solutions = false;

        /**
        * The repair function applied to each Candidate of the population after the mutations if it isn't a nullptr. \n
        * This can be used to perform local search after the mutations, implementing a memetic algorithm.
        */
        repairFunction_t repairFunction = nullptr;

        callbackFunction_t endOfGenerationCallback = nullptr;

        /**
        * Standard constructor for the GA.
        *
        * @param chrom_len The length of the chromosomes (number of genes).
        * @param fitness_function The fitness function to find the maximum of.
        */
        GA(size_t chrom_len, fitnessFunction_t fitness_function);

        virtual ~GA() = default;

        /**
        * Runs the genetic algorithm with the selected settings.
        *
        * @returns The optimal solutions.
        */
        [[maybe_unused]] CandidateVec run();


        /** @returns A vector of the pareto optimal solutions found while running the algorithm. */
        [[nodiscard]] CandidateVec solutions() const;

        /** @returns The number of fitness evaluations performed while running the algorithm. */
        [[nodiscard]] size_t num_fitness_evals() const;

        /** @returns The current value of the generation counter. */
        [[nodiscard]] size_t generation_cntr() const;

        /** @returns The population of the final generation in the algorithm. */
        [[nodiscard]] Population population() const;

        /** @returns A History object containing stats from each generation of the single objective genetic algorithm. */
        [[nodiscard]] History soga_history() const;

        /**
        * Set the type of the problem/genetic algorithm that will be used (single-/multi-objective).
        *
        * @param mode The algorithm type to use. @see Mode
        */
        void mode(Mode mode);
        [[nodiscard]] Mode mode() const;

        /**
        * Sets the length of the chromosomes (number of genes) of the Candidate solutions used in the algorithm to @p len. \n
        * The chromosome length must be at least 1.
        *
        * @param len The length of the chromosomes.
        */
        void chrom_len(size_t len);
        [[nodiscard]] size_t chrom_len() const;

        /**
        * Sets the number of Candidates used in the population to @p size. \n
        * The population size must be at least 1.
        *
        * @param size The size of the populations.
        */
        void population_size(size_t size);
        [[nodiscard]] size_t population_size() const;

        /**
        * Sets the crossover rate used in the algorithm to @p pc. \n
        * The crossover rate must be on the closed interval [0.0, 1.0].
        *
        * @param pc The probability of crossover being performed on 2 selected individuals.
        */
        void crossover_rate(double pc);
        [[nodiscard]] double crossover_rate() const;

        /**
        * Sets the mutation rate used in the algorithm to @p pm. \n
        * The mutation rate must be on the closed interval [0.0, 1.0].
        *
        * @param pm The probability of mutating the children generated by the crossovers (per child or per gene based on encoding type).
        */
        void mutation_rate(double pm);
        [[nodiscard]] double mutation_rate() const;

        /**
        * Sets the selection function used in the single-objective algorithm to @p f. \n
        * The selection function set is ignored in the other algorithm types. @see Mode
        *
        * @param f The selection function used in the single-objective algorithm.
        */
        void selection_method(selectionFunction_t f);

        /**
        * Sets the selection method used in the single-objective algorithm to @p method. \n
        * The selection method set is ignored in the other algorithm types. @see Mode
        *
        * @param method The selection method used in the single-objective algorithm.
        */
        void selection_method(SogaSelection method);
        [[nodiscard]] SogaSelection selection_method() const;

        /**
        * Sets the number of solutions picked for the tournaments to @p size if the tournament selection method is
        * selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The size of the tournaments must be at least 2.
        *
        * @param size The size of the tournaments during tournament selection.
        */
        void tournament_size(size_t size);
        [[nodiscard]] size_t tournament_size() const;

        /**
        * Sets the minimum and maximum selection weights used during rank selection if the rank selection method is
        * selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The @p min_weight must be on the closed interval [0.0, @p max_weight]. \n
        * The @p max_weight must be greater than @p min_weight.
        *
        * @param min_weight The selection weight assigned to the worst Candidate in the population.
        * @param max_weight The selection weight assigned to the best Candidate in the population.
        */
        void rank_sel_weights(double min_weight, double max_weight);
        [[nodiscard]] std::pair<double, double> rank_sel_weights() const;

        /**
        * Sets the minimum and maximum temperature values used during boltzmann selection if the boltzmann
        * selection method is selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The minimum @p tmin temperature must be on the interval [0.1, @p tmax). \n
        * The maximum @p tmax temperature must be on the interval (@p tmin, DBL_MAX).
        *
        * @param tmin The minimum temperature (used in the last generation).
        * @param tmax The maximum temperature (used in the first generation).
        */
        void boltzmann_temps(double tmin, double tmax);
        [[nodiscard]] std::pair<double, double> boltzmann_temps() const;

        /**
        * Sets the scaling parameter used during selection to @p scale if the sigma selection method is
        * selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The scaled fitness values are: f' = (f - f_mean) / (@p scale * f_sd) \n
        * The value of @p scale must be on the interval [1.0, DBL_MAX].
        *
        * @param scale The scaling parameter used during sigma selection.
        */
        void sigma_scale(double scale);
        [[nodiscard]] double sigma_scale() const;

        /**
        * Sets the stop condition used in the algorithm to @p condition. Some of the stop conditions
        * only work with the single-objective algorithm. @see StopCondition \n
        * The algorithm always stops when the set max_gen generation has been reached regardless of the stop condition
        * set. @see max_gen
        *
        * @param condition The stop condition used in the algorithm.
        */
        void stop_condition(StopCondition condition);
        [[nodiscard]] StopCondition stop_condition() const;

        /**
        * Sets the maximum number of generations the algorithm runs for to @p max_gen. The
        * algorithm will always stop when this generation has been reached regardless of what stop
        * condition was set, but it can stop earlier when another stop condition is selected.
        * @see stop_condition @see StopCondition \n
        * The value of @p max_gen must be at least 1.
        *
        * @param max_gen The maximum number of generations.
        */
        void max_gen(size_t max_gen);
        [[nodiscard]] size_t max_gen() const;

        /**
        * Sets the maximum number of fitness evaluations the algorithm runs for to @p max_evals if
        * the fitness_evals stop condition is selected. @see stop_condition @see StopCondition \n
        * The algorithm may evaluate the fitness function more times than the maximum set since the
        * stop condition is only checked at the end of each generation. \n
        * The value of @p max_evals must be at least 1.
        *
        * @param max_evals
        */
        void max_fitness_evals(size_t max_evals);
        [[nodiscard]] size_t max_fitness_evals() const;

        /**
        * Sets the reference fitness value for the fitness_value stop condition to @p ref. \n
        * The algorithm will stop running if a solution has been found which dominates this reference point. \n
        * The size of the reference vector should be equal to the number of objectives. \n
        * @see stop_condition @see StopCondition
        *
        * @param ref The fitness reference used.
        */
        void fitness_threshold(std::vector<double> ref);
        [[nodiscard]] std::vector<double> fitness_threshold() const;

        /**
        * Sets the number of generations to look back when evaluating the stall stop conditions. \n
        * Only relevant for the single-objective algorithm. @see stop_condition @see StopCondition \n
        * Must be at least 1.
        *
        * @param count The number of generations to look back when checking the stall conditions.
        */
        void stall_gen_count(size_t count);
        [[nodiscard]] size_t stall_gen_count() const;

        /**
        * Sets the value of the stall threshold to @p threshold for the stall stop conditions. \n
        * Only relevant for the single-objective algorithm. @see stop_condition @see StopCondition \n
        * May be negative if the deterioration of the stall metric can be allowed.
        *
        * @param threshold The stall threshold to use.
        */
        void stall_threshold(double threshold);
        [[nodiscard]] double stall_threshold() const;

        /**
        * Sets the initial population to be used in the algorithm to @p pop instead of randomly generating it. \n
        * If @p pop is empty, the initial population will be randomly generated. \n
        * If the preset population's size is not equal to the population size set, either additional randomly generated
        * Candidates will be added to fill out the initial population, or some Candidates will be discarded from the preset.
        *
        * @param pop The initial population to use in the algorithm.
        */
        void presetInitialPopulation(const Population& pop);

        /**
        * Sets the fitness function used by the algorithm to @p f. \n
        * The fitness function should return a vector whose size is equal to the number of objectives, and
        * each element of the vector should be finite.
        *
        * @param f The fitness function to find the maximum of.
        */
        void setFitnessFunction(fitnessFunction_t f);

        /* Some getters for the NSGA-III algorithm. */
        [[nodiscard]] std::vector<std::vector<double>> ref_points() const;
        [[nodiscard]] std::vector<double> ideal_point() const;
        [[nodiscard]] std::vector<double> nadir_point() const;

    protected:

        Population population_;
        size_t generation_cntr_ = 0;
        size_t num_objectives_ = 0;		/* Determined from the fitness function. */

        /* For the NSGA-III. */
        std::vector<std::vector<double>> ref_points_;
        std::vector<double> ideal_point_;
        std::vector<double> nadir_point_;
        std::vector<std::vector<double>> extreme_points_;

        /* Results of the GA. */
        CandidateVec solutions_;
        std::atomic<size_t> num_fitness_evals_ = 0;
        History soga_history_;

        /* Basic parameters of the GA. */
        Mode mode_ = Mode::single_objective;
        size_t chrom_len_;
        size_t population_size_ = 100;
        double crossover_rate_ = 0.8;
        double mutation_rate_ = 0.01;

        /* Single-objective GA selection settings. */
        SogaSelection selection_method_ = SogaSelection::tournament;
        size_t tournament_size_ = 2;
        double rank_sel_min_w_ = 0.1;
        double rank_sel_max_w_ = 1.1;
        double boltzmann_tmin_ = 0.25;
        double boltzmann_tmax_ = 4.0;
        double sigma_scale_ = 3.0;

        /* Stop condition settings. */
        StopCondition stop_condition_ = StopCondition::max_gen;
        size_t max_gen_ = 500;
        size_t max_fitness_evals_ = 5000;
        std::vector<double> fitness_reference_;
        size_t stall_gen_count_ = 20;
        double stall_threshold_ = 1e-6;

        /* Initial population settings. */
        Population initial_population_preset_;

        /* User supplied functions used in the GA. All of these are optional except for the fitness function. */
        fitnessFunction_t fitnessFunction;
        selectionFunction_t customSelection = nullptr;
        crossoverFunction_t customCrossover = nullptr;
        mutationFunction_t customMutate = nullptr;


        /* General functions for the genetic algorithms. */

        void init();
        virtual Candidate generateCandidate() const = 0;
        Population generateInitialPopulation() const;
        void evaluate(Population& pop);
        void updateOptimalSolutions(CandidateVec& optimal_sols, const Population& pop) const;
        void prepSelections(Population& pop) const;
        Candidate select(const Population& pop) const;
        virtual CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const = 0;
        virtual void mutate(Candidate& child) const = 0;
        void repair(Population& pop) const;      
        Population updatePopulation(Population& old_pop, CandidateVec& children);       
        bool stopCondition() const;
        void updateStats(const Population& pop);


        /* SOGA functions. */

        /* Functions for calculating the selection probabilities of individuals in the single-objective algorithm. */

        static void sogaCalcRouletteWeights(Population& pop);
        static void sogaCalcRankWeights(Population& pop, double weight_min = 0.1, double weight_max = 1.1);
        static void sogaCalcSigmaWeights(Population& pop, double scale = 3.0);
        static void sogaCalcBoltzmannWeights(Population& pop, size_t t, size_t t_max, double temp_min, double temp_max);

        void sogaCalcWeights(Population& pop) const;

        /* Functions used for the selections in the single-objective algorithm. */

        static Candidate sogaWeightProportionalSelect(const Population& pop);
        static Candidate sogaTournamentSelect(const Population& pop, size_t tourney_size);

        Candidate sogaSelect(const Population& pop) const;

        /* Create the population of the next generation from the old population and the children. */
        Population updateSogaPopulation(Population& old_pop, CandidateVec& children) const;


        /* NSGA-II functions. */

        /* Find all Pareto fronts in the population and also assign the nondomination ranks of the candidates (assuming fitness maximization). */
        static std::vector<std::vector<size_t>> nonDominatedSort(Population& pop);

        /* Calculate the crowding distances of the candidates in each pareto front in pfronts of the population. */
        static void calcCrowdingDistances(Population& pop, std::vector<std::vector<size_t>>& pfronts);

        /* Returns true if lhs is better than rhs. */
        static bool crowdedCompare(const Candidate& lhs, const Candidate& rhs);

        static Candidate nsga2Select(const Population& pop);

        /* Create the population of the next generation from the old population and the children. */
        Population updateNsga2Population(Population& old_pop, CandidateVec& children) const;


        /* NSGA-III functions. */

        void updateIdealPoint(const Population& pop);
        void updateNadirPoint(const Population& pop);

        /* Find the closest reference point to each candidate after normalization, and their distances. */
        void associatePopToRefs(Population& pop, const std::vector<std::vector<double>>& ref_points);

        /* Return the niche counts of the ref points and assign niche counts to the candidates. */
        static std::vector<size_t> calcNicheCounts(Population& pop, const std::vector<std::vector<double>>& ref_points);

        /* Returns true if lhs is better than rhs. */
        static bool nichedCompare(const Candidate& lhs, const Candidate& rhs);

        /* Tournament selection using the niche counts for tiebreaks. */
        static Candidate nsga3Select(const Population& pop);

        /* Create the population of the next generation from the old population and the children. */
        Population updateNsga3Population(Population& old_pop, CandidateVec& children);

        /* Functions used to find the optimal solutions in the population. */
        static CandidateVec findParetoFront1D(const Population& pop);
        static CandidateVec findParetoFrontKung(const Population& pop);


        /* Utility functions. */

        /* Find the minimum/maximum fitness along each objective in the population. */
        static std::vector<double> fitnessMin(const Population& pop);
        static std::vector<double> fitnessMax(const Population& pop);

        /* Find the mean/standard deviation of the fitness values of the population along the first objective. */
        static double fitnessMean(const Population& pop);
        static double fitnessSD(const Population& pop);

    };

} // namespace genetic_algorithm


/* IMPLEMENTATION */

#include <execution>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <cassert>
#include <cmath>

#include "rng.h"
#include "reference_points.h"
#include "mo_detail.h"

namespace genetic_algorithm
{
    template<typename geneType>
    inline void GA<geneType>::History::clear() noexcept
    {
        fitness_mean.clear();
        fitness_sd.clear();
        fitness_min.clear();
        fitness_max.clear();
    }

    template<typename geneType>
    inline void GA<geneType>::History::reserve(size_t new_capacity)
    {
        fitness_mean.reserve(new_capacity);
        fitness_sd.reserve(new_capacity);
        fitness_min.reserve(new_capacity);
        fitness_max.reserve(new_capacity);
    }

    template<typename geneType>
    inline void GA<geneType>::History::add(double mean, double sd, double min, double max)
    {
        fitness_mean.push_back(mean);
        fitness_sd.push_back(sd);
        fitness_min.push_back(min);
        fitness_max.push_back(max);
    }

    template<typename geneType>
    inline GA<geneType>::Candidate::Candidate()
    {
    }

    template<typename geneType>
    inline GA<geneType>::Candidate::Candidate(const std::vector<geneType>& chrom)
        : chromosome(chrom)
    {
    }

    template<typename geneType>
    inline GA<geneType>::Candidate::Candidate(std::vector<geneType>&& chrom) noexcept
        : chromosome(std::move(chrom))
    {
    }

    template<typename geneType>
    inline bool GA<geneType>::Candidate::operator==(const Candidate& rhs) const
    {
        return (this->chromosome == rhs.chromosome);
    }

    template<typename geneType>
    inline bool GA<geneType>::Candidate::operator!=(const Candidate& rhs) const
    {
        return !(this->chromosome == rhs.chromosome);
    }

    template<typename geneType>
    inline size_t GA<geneType>::CandidateHasher::operator()(const Candidate& c) const noexcept
    {
        size_t seed = c.chromosome.size();
        for (const auto& gene : c.chromosome)
        {
            seed ^= std::hash<geneType>()(gene) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    template<typename geneType>
    inline GA<geneType>::GA(size_t chrom_len, fitnessFunction_t fitness_function)
        : chrom_len_(chrom_len), mutation_rate_(1.0 / chrom_len), fitnessFunction(fitness_function)
    {
        if (chrom_len == 0)
        {
            throw std::invalid_argument("The chromosome length must be at least 1.");
        }
        if (fitnessFunction == nullptr)
        {
            throw std::invalid_argument("The fitness function is a nullptr.");
        }
    }

    template<typename geneType>
    inline typename GA<geneType>::CandidateVec GA<geneType>::solutions() const
    {
        return solutions_;
    }

    template<typename geneType>
    inline size_t GA<geneType>::num_fitness_evals() const
    {
        return static_cast<size_t>(num_fitness_evals_);
    }

    template<typename geneType>
    inline size_t GA<geneType>::generation_cntr() const
    {
        return generation_cntr_;
    }

    template<typename geneType>
    inline typename GA<geneType>::Population GA<geneType>::population() const
    {
        return population_;
    }

    template<typename geneType>
    inline typename GA<geneType>::History GA<geneType>::soga_history() const
    {
        return soga_history_;
    }

    template<typename geneType>
    inline void GA<geneType>::mode(Mode mode)
    {
        if (static_cast<size_t>(mode) > 2) throw std::invalid_argument("Invalid algorithm mode selected.");

        mode_ = mode;
    }

    template<typename geneType>
    inline typename GA<geneType>::Mode GA<geneType>::mode() const
    {
        return mode_;
    }

    template<typename geneType>
    inline void GA<geneType>::chrom_len(size_t len)
    {
        if (len == 0) throw std::invalid_argument("The chromosome length must be at least 1.");

        chrom_len_ = len;
    }

    template<typename geneType>
    inline size_t GA<geneType>::chrom_len() const
    {
        return chrom_len_;
    }

    template<typename geneType>
    inline void GA<geneType>::population_size(size_t size)
    {
        if (size == 0) throw std::invalid_argument("The population size must be at least 1.");

        population_size_ = size;
    }

    template<typename geneType>
    inline size_t GA<geneType>::population_size() const
    {
        return population_size_;
    }

    template<typename geneType>
    inline void GA<geneType>::crossover_rate(double pc)
    {
        if (!(0.0 <= pc && pc <= 1.0)) throw std::invalid_argument("The crossover probability must be in the range [0, 1].");

        crossover_rate_ = pc;
    }

    template<typename geneType>
    inline double GA<geneType>::crossover_rate() const
    {
        return crossover_rate_;
    }

    template<typename geneType>
    inline void GA<geneType>::mutation_rate(double pm)
    {
        if (!(0.0 <= pm && pm <= 1.0)) throw std::invalid_argument("The mutation probability must be in the range [0, 1].");

        mutation_rate_ = pm;
    }

    template<typename geneType>
    inline double GA<geneType>::mutation_rate() const
    {
        return mutation_rate_;
    }

    template<typename geneType>
    inline void GA<geneType>::selection_method(selectionFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The selection function can't be a nullptr.");

        selection_method_ = SogaSelection::custom;
        customSelection = f;
    }

    template<typename geneType>
    inline void GA<geneType>::selection_method(SogaSelection method)
    {
        if (static_cast<size_t>(method) > 5) throw std::invalid_argument("Invalid soga selection method selected.");

        selection_method_ = method;
    }

    template<typename geneType>
    inline typename GA<geneType>::SogaSelection GA<geneType>::selection_method() const
    {
        return selection_method_;
    }

    template<typename geneType>
    inline void GA<geneType>::tournament_size(size_t size)
    {
        if (size < 2) throw std::invalid_argument("The tournament size must be at least 2.");

        tournament_size_ = size;
    }

    template<typename geneType>
    inline size_t GA<geneType>::tournament_size() const
    {
        return tournament_size_;
    }

    template<typename geneType>
    inline void GA<geneType>::rank_sel_weights(double min_weight, double max_weight)
    {
        if (!(0.0 <= min_weight && min_weight <= max_weight))
        {
            throw std::invalid_argument("The minimum weight must be in the range [0.0, max_weight].");
        }
        if (!(min_weight <= max_weight && max_weight <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("The maximum weight must be in the range [min_weight, DBL_MAX].");
        }

        rank_sel_min_w_ = min_weight;
        rank_sel_max_w_ = max_weight;
    }

    template<typename geneType>
    inline std::pair<double, double> GA<geneType>::rank_sel_weights() const
    {
        return { rank_sel_min_w_, rank_sel_max_w_ };
    }

    template<typename geneType>
    inline void GA<geneType>::boltzmann_temps(double tmin, double tmax)
    {
        if (!(0.1 <= tmin && tmin < tmax))
        {
            throw std::invalid_argument("The minimum temperature (tmin) must be in the range [0.1, tmax).");
        }
        if (!(tmin < tmax && tmax <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("The maximum temperature (tmax) must be in the range (tmin, DBL_MAX].");
        }

        boltzmann_tmin_ = tmin;
        boltzmann_tmax_ = tmax;
    }

    template<typename geneType>
    inline std::pair<double, double> GA<geneType>::boltzmann_temps() const
    {
        return { boltzmann_tmin_, boltzmann_tmax_ };
    }

    template<typename geneType>
    inline void GA<geneType>::sigma_scale(double scale)
    {
        if (!(1.0 <= scale && scale <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Scale must be in the range [1.0, DBL_MAX].");
        }

        sigma_scale_ = scale;
    }

    template<typename geneType>
    inline double GA<geneType>::sigma_scale() const
    {
        return sigma_scale_;
    }

    template<typename geneType>
    inline void GA<geneType>::stop_condition(StopCondition condition)
    {
        if (static_cast<size_t>(condition) > 4) throw std::invalid_argument("Invalid stop condition selected.");

        stop_condition_ = condition;
    }

    template<typename geneType>
    inline typename GA<geneType>::StopCondition GA<geneType>::stop_condition() const
    {
        return stop_condition_;
    }

    template<typename geneType>
    inline void GA<geneType>::max_gen(size_t max_gen)
    {
        if (max_gen == 0) throw std::invalid_argument("The maximum number of generations must be at least 1.");

        max_gen_ = max_gen;
    }

    template<typename geneType>
    inline size_t GA<geneType>::max_gen() const
    {
        return max_gen_;
    }

    template<typename geneType>
    inline void GA<geneType>::max_fitness_evals(size_t max_evals)
    {
        if (max_evals == 0) throw std::invalid_argument("The maximum number of fitness evaluations must be at least 1.");

        max_fitness_evals_ = max_evals;
    }

    template<typename geneType>
    inline size_t GA<geneType>::max_fitness_evals() const
    {
        return max_fitness_evals_;
    }

    template<typename geneType>
    inline void GA<geneType>::fitness_threshold(std::vector<double> ref)
    {
        if (ref.empty())
        {
            throw std::invalid_argument("The reference vector is empty.");
        }
        if (!std::all_of(ref.begin(), ref.end(), [](double val) { return std::isfinite(val); }))
        {
            throw std::invalid_argument("Invalid value in the reference vector.");
        }

        fitness_reference_ = ref;
    }

    template<typename geneType>
    inline std::vector<double> GA<geneType>::fitness_threshold() const
    {
        return fitness_reference_;
    }

    template<typename geneType>
    inline void GA<geneType>::stall_gen_count(size_t count)
    {
        if (count == 0) throw std::invalid_argument("The stall generation count must be at least 1.");

        stall_gen_count_ = count;
    }

    template<typename geneType>
    inline size_t GA<geneType>::stall_gen_count() const
    {
        return stall_gen_count_;
    }

    template<typename geneType>
    inline void GA<geneType>::stall_threshold(double threshold)
    {
        if (!std::isfinite(threshold)) throw std::invalid_argument("The stall threshold must be finite.");

        stall_threshold_ = threshold;
    }

    template<typename geneType>
    inline double GA<geneType>::stall_threshold() const
    {
        return stall_threshold_;
    }

    template<typename geneType>
    inline void GA<geneType>::presetInitialPopulation(const Population& pop)
    {
        if (!std::all_of(pop.begin(), pop.end(), [this](const Candidate& c) { return c.chromosome.size() == chrom_len_; }))
        {
            throw std::invalid_argument("The length of each chromosome in the preset pop must be equal to chrom_len.");
        }

        initial_population_preset_ = pop;
    }

    template<typename geneType>
    inline void GA<geneType>::setFitnessFunction(fitnessFunction_t f)
    {
        if (f == nullptr) throw std::invalid_argument("The fitness function can't be a nullptr.");

        fitnessFunction = f;
    }

    template<typename geneType>
    inline std::vector<std::vector<double>> GA<geneType>::ref_points() const
    {
        return ref_points_;
    }

    template<typename geneType>
    inline std::vector<double> GA<geneType>::ideal_point() const
    {
        return ideal_point_;
    }

    template<typename geneType>
    inline std::vector<double> GA<geneType>::nadir_point() const
    {
        return nadir_point_;
    }


    template<typename geneType>
    inline typename GA<geneType>::CandidateVec GA<geneType>::run()
    {
        using namespace std;

        init();

        /* Create and evaluate the initial population. */
        population_ = generateInitialPopulation();
        evaluate(population_);
        updateStats(population_);

        /* Other generations. */
        size_t num_children = population_size_ + population_size_ % 2;
        while (!stopCondition())
        {
            vector<CandidatePair> parent_pairs(num_children / 2);

            prepSelections(population_);
            if (archive_optimal_solutions) updateOptimalSolutions(solutions_, population_);

            /* Selections. */
            generate(execution::par_unseq, parent_pairs.begin(), parent_pairs.end(),
            [this]() -> CandidatePair
            {
                return make_pair(select(population_), select(population_));
            });

            /* Crossovers. */
            for_each(execution::par_unseq, parent_pairs.begin(), parent_pairs.end(),
            [this](CandidatePair& p) -> void
            {
                p = crossover(p.first, p.second);
            });

            vector<Candidate> children;
            children.reserve(num_children);
            for (size_t i = 0; i < parent_pairs.size(); i++)
            {
                children.push_back(move(parent_pairs[i].first));
                children.push_back(move(parent_pairs[i].second));
            }

            /* Mutations. */
            for_each(execution::par_unseq, children.begin(), children.end(),
            [this](Candidate& c) -> void
            {
                mutate(c);
            });

            /* Apply repair function to the children if set. */
            repair(children);

            /* Overwrite the current population with the children. */
            evaluate(children);
            population_ = updatePopulation(population_, children);

            if (endOfGenerationCallback != nullptr) endOfGenerationCallback(this);
            generation_cntr_++;

            updateStats(population_);
        }
        updateOptimalSolutions(solutions_, population_);

        return solutions_;
    }

    template<typename geneType>
    inline void GA<geneType>::init()
    {
        /* Check stop condition. */
        if (mode_ != Mode::single_objective)
        {
            if (stop_condition_ == StopCondition::fitness_mean_stall || stop_condition_ == StopCondition::fitness_best_stall)
            {
                throw std::invalid_argument("The stall stop conditions only work for the single-objective algorithm.");
            }
        }
        /* Check selection method. */
        if (selection_method_ == SogaSelection::custom && customSelection == nullptr)
        {
            throw std::invalid_argument("The custom selection function is a nullptr.");
        }
        /* Check mode. */
        Candidate temp = generateCandidate();
        temp.fitness = fitnessFunction(temp.chromosome);
        num_objectives_ = temp.fitness.size();
        if (mode_ == Mode::single_objective && num_objectives_ != 1)
        {
            throw std::invalid_argument("The size of the fitness vector must be 1 for single-objective optimization.");
        }
        else if (mode_ != Mode::single_objective && num_objectives_ < 2)
        {
            throw std::invalid_argument("The size of the fitness vector must be at least 2 for multi-objective optimization.");
        }

        /* General initialization. */
        generation_cntr_ = 0;
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();

        /* Single objective stuff. */
        if (mode_ == Mode::single_objective)
        {
            soga_history_.clear();
            soga_history_.reserve(max_gen_);
        }

        /* Multi-objective stuff (NSGA-III). */
        ideal_point_ = std::vector<double>(num_objectives_, -std::numeric_limits<double>::max());
        nadir_point_ = std::vector<double>(num_objectives_);
        extreme_points_ = std::vector<std::vector<double>>(num_objectives_, std::vector<double>(num_objectives_));

        /* Generate the reference points for the NSGA-III algorithm. */
        if (mode_ == Mode::multi_objective_decomp)
        {
            ref_points_ = detail::generateRefPoints(population_size_, num_objectives_);
        }
    }

    template<typename geneType>
    inline typename GA<geneType>::Population GA<geneType>::generateInitialPopulation() const
    {
        assert(population_size_ > 0);

        if (!std::all_of(initial_population_preset_.begin(), initial_population_preset_.end(),
        [this](const Candidate& sol)
        {
            return sol.chromosome.size() == chrom_len_;
        }))
        {
            throw std::length_error("The chromosome lengths in the preset initial population must be equal to the chrom_len set.");
        }

        Population pop;
        pop.reserve(population_size_);

        for (size_t i = 0; i < std::min(population_size_, initial_population_preset_.size()); i++)
        {
            pop.push_back(initial_population_preset_[i]);
        }
        while (pop.size() < population_size_)
        {
            pop.push_back(generateCandidate());
        }

        return pop;
    }

    template<typename geneType>
    inline void GA<geneType>::evaluate(Population& pop)
    {
        assert(fitnessFunction != nullptr);

        std::for_each(std::execution::par_unseq, pop.begin(), pop.end(),
        [this](Candidate& sol)
        {
            if (changing_fitness_func || !sol.is_evaluated)
            {
                sol.fitness = fitnessFunction(sol.chromosome);
                sol.is_evaluated = true;

                num_fitness_evals_++;
            }
        });

        for (const auto& sol : pop)
        {
            if (sol.fitness.size() != num_objectives_)
            {
                throw std::domain_error("A fitness vector returned by the fitness function has incorrect size.");
            }
            if (!std::all_of(sol.fitness.begin(), sol.fitness.end(), [](double val) { return std::isfinite(val); }))
            {
                throw std::domain_error("A non-finite fitness value was returned by the fitness function.");
            }
        }
    }

    template<typename geneType>
    inline void GA<geneType>::updateOptimalSolutions(CandidateVec& optimal_sols, const Population& pop) const
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.is_evaluated; }));

        optimal_sols.insert(optimal_sols.end(), pop.begin(), pop.end());
        if (mode_ == Mode::single_objective)
        {
            optimal_sols = findParetoFront1D(optimal_sols);
        }
        else
        {
            optimal_sols = findParetoFrontKung(optimal_sols);
        }

        /* Remove duplicate solutions. */
        CandidateSet unique_sols;
        for (const auto& sol : optimal_sols) unique_sols.insert(std::move(sol));	/* Insertion is faster than ctor. */
        optimal_sols.assign(unique_sols.begin(), unique_sols.end());
    }

    template<typename geneType>
    inline void GA<geneType>::prepSelections(Population& pop) const
    {
        switch (mode_)
        {
            case Mode::single_objective:
                sogaCalcWeights(pop);
                break;
            case Mode::multi_objective_sorting:
                /* Nothing to do. */
                break;
            case Mode::multi_objective_decomp:
                /* Nothing to do. */
                break;
            default:
                assert(false);	/* Invalid mode. Shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline typename GA<geneType>::Candidate GA<geneType>::select(const Population& pop) const
    {
        switch (mode_)
        {
            case Mode::single_objective:
                return sogaSelect(pop);
            case Mode::multi_objective_sorting:
                return nsga2Select(pop);
            case Mode::multi_objective_decomp:
                return nsga3Select(pop);
            default:
                assert(false);	/* Invalid mode. Shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline void GA<geneType>::repair(Population& pop) const
    {
        /* This function doesn't do anything unless a repair function is specified. */
        if (repairFunction == nullptr) return;

        std::for_each(std::execution::par_unseq, pop.begin(), pop.end(),
        [this](Candidate& sol)
        {
            Chromosome improved_chrom = repairFunction(sol.chromosome);
            if (improved_chrom != sol.chromosome)
            {
                sol.is_evaluated = false;
                sol.chromosome = std::move(improved_chrom);
            }
        });

        for (const auto& sol : pop)
        {
            if (sol.chromosome.size() != chrom_len_)
            {
                throw std::domain_error("The repair function must return chromosomes of chrom_len length.");
            }
        }
    }

    template<typename geneType>
    inline typename GA<geneType>::Population GA<geneType>::updatePopulation(Population& old_pop, CandidateVec& children)
    {
        switch (mode_)
        {
            case Mode::single_objective:
                return updateSogaPopulation(old_pop, children);
            case Mode::multi_objective_sorting:
                return updateNsga2Population(old_pop, children);
            case Mode::multi_objective_decomp:
                return updateNsga3Population(old_pop, children);
            default:
                assert(false);	/* Invalid mode, shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline bool GA<geneType>::stopCondition() const
    {
        if (mode_ != Mode::single_objective && stop_condition_ == StopCondition::fitness_best_stall)
        {
            throw std::invalid_argument("The stall stop conditions only work with the single-objective algorithm.");
        }
        else if (mode_ != Mode::single_objective && stop_condition_ == StopCondition::fitness_mean_stall)
        {
            throw std::invalid_argument("The stall stop conditions only work with the single-objective algorithm.");
        }

        /* Always stop when reaching max_gen regardless of stop condition. */
        if (generation_cntr_ >= max_gen_ - 1) return true;

        /* Early-stop conditions. */
        double metric_now, metric_old;
        switch (stop_condition_)
        {
            case StopCondition::max_gen:
                /* Already checked above. */
                return false;

            case StopCondition::fitness_value:
                return std::any_of(population_.begin(), population_.end(),
                [this](const Candidate& sol)
                {
                    return detail::paretoCompare(fitness_reference_, sol.fitness);
                });

            case StopCondition::fitness_evals:
                return num_fitness_evals_ >= max_fitness_evals_;

            case StopCondition::fitness_mean_stall:
                if (generation_cntr_ >= stall_gen_count_)
                {
                    metric_now = soga_history_.fitness_mean[generation_cntr_];
                    metric_old = soga_history_.fitness_mean[generation_cntr_ - stall_gen_count_];

                    return (metric_now - metric_old) < stall_threshold_;
                }
                else return false;

            case StopCondition::fitness_best_stall:
                if (generation_cntr_ >= stall_gen_count_)
                {
                    metric_now = soga_history_.fitness_mean[generation_cntr_];
                    metric_old = soga_history_.fitness_mean[generation_cntr_ - stall_gen_count_];

                    return (metric_now - metric_old) < stall_threshold_;
                }
                else return false;

            default:
                assert(false);	/* Invalid stop condition. Shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline void GA<geneType>::updateStats(const Population& pop)
    {
        switch (mode_)
        {
            case Mode::single_objective:
                soga_history_.add(fitnessMean(pop), fitnessSD(pop), fitnessMin(pop)[0], fitnessMax(pop)[0]);
                break;
            case Mode::multi_objective_sorting:
                break;
            case Mode::multi_objective_decomp:
                break;
            default:
                assert(false);	/* Invalid mode, shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline void GA<geneType>::sogaCalcRouletteWeights(Population& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        /* Roulette selection wouldn't work for negative fitness values. */
        bool has_negative_fitness = std::any_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness[0] < 0.0; });
        double offset = fitnessMin(pop)[0] * has_negative_fitness;

        double pdf_mean = 0.0;
        for (auto& sol : pop)
        {
            sol.selection_pdf = sol.fitness[0] - 2.0 * offset;
            pdf_mean += sol.selection_pdf / pop.size();
        }

        double cdf = 0.0;
        for (auto& sol : pop)
        {
            sol.selection_pdf /= pdf_mean;
            sol.selection_pdf /= pop.size();

            cdf += sol.selection_pdf;
            sol.selection_cdf = cdf;
        }
    }

    template<typename geneType>
    inline void GA<geneType>::sogaCalcRankWeights(Population& pop, double weight_min, double weight_max)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
        assert(0.0 <= weight_min && weight_min < weight_max&& weight_max <= std::numeric_limits<double>::max());

        /* Argsort descending order. */
        std::vector<size_t> indices(pop.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(),
        [&pop](size_t lidx, size_t ridx)
        {
            return pop[lidx].fitness[0] > pop[ridx].fitness[0];
        });

        double pdf_mean = 0.0;
        for (size_t i = 0; i < indices.size(); i++)
        {
            double m = 1.0 - i / (pop.size() - 1.0);
            pop[indices[i]].selection_pdf = weight_min + (weight_max - weight_min) * m;

            pdf_mean += pop[indices[i]].selection_pdf / pop.size();
        }

        double cdf = 0.0;
        for (auto& sol : pop)
        {
            sol.selection_pdf /= pdf_mean;
            sol.selection_pdf /= pop.size();

            cdf += sol.selection_pdf;
            sol.selection_cdf = cdf;
        }
    }

    template<typename geneType>
    inline void GA<geneType>::sogaCalcSigmaWeights(Population& pop, double scale)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
        assert(scale > 1.0);

        double fitness_mean = fitnessMean(pop);
        double fitness_sd = fitnessSD(pop);

        double pdf_mean = 0.0;
        for (auto& sol : pop)
        {
            sol.selection_pdf = 1.0 + (sol.fitness[0] - fitness_mean) / (scale * std::max(fitness_sd, 1E-6));
            
            /* If (fitness < f_mean - scale * SD) the weight could be negative. */
            sol.selection_pdf = std::max(sol.selection_pdf, 0.0);

            pdf_mean += sol.selection_pdf / pop.size();
        }

        double cdf = 0.0;
        for (auto& sol : pop)
        {
            sol.selection_pdf /= pdf_mean;
            sol.selection_pdf /= pop.size();

            cdf += sol.selection_pdf;
            sol.selection_cdf = cdf;
        }
    }

    template<typename geneType>
    inline void GA<geneType>::sogaCalcBoltzmannWeights(Population& pop, size_t t, size_t t_max, double temp_min, double temp_max)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
        assert(t_max >= t);
        assert(temp_max > temp_min && temp_min > 0.1);

        double temperature = -temp_max / (1.0 + std::exp(-10.0 * (double(t) / t_max) + 3.0)) + temp_max + temp_min;

        double fmax = fitnessMax(pop)[0];
        double fmin = fitnessMin(pop)[0];

        double pdf_mean = 0.0;
        for (auto& sol : pop)
        {
            /* Norm fitness values so the exp function won't return too high values. */
            double fnorm = (sol.fitness[0] - fmin) / std::max(fmax - fmin, 1E-6);

            sol.selection_pdf = exp(fnorm / temperature);
            pdf_mean += sol.selection_pdf / pop.size();
        }

        double cdf = 0.0;
        for (auto& sol : pop)
        {
            sol.selection_pdf /= pdf_mean;
            sol.selection_pdf /= pop.size();

            cdf += sol.selection_pdf;
            sol.selection_cdf = cdf;
        }
    }

    template<typename geneType>
    inline void GA<geneType>::sogaCalcWeights(Population& pop) const
    {
        switch (selection_method_)
        {
            case SogaSelection::tournament:
                /* Not needed for tournament selection. */
                break;
            case SogaSelection::roulette:
                sogaCalcRouletteWeights(pop);
                break;
            case SogaSelection::rank:
                sogaCalcRankWeights(pop, rank_sel_min_w_, rank_sel_max_w_);
                break;
            case SogaSelection::sigma:
                sogaCalcSigmaWeights(pop, sigma_scale_);
                break;
            case SogaSelection::boltzmann:
                sogaCalcBoltzmannWeights(pop, generation_cntr_, max_gen_, boltzmann_tmin_, boltzmann_tmax_);
                break;
            case SogaSelection::custom:
                break;
            default:
                assert(false);	/* Invalid selection method. Shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline typename GA<geneType>::Candidate GA<geneType>::sogaWeightProportionalSelect(const Population& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        double threshold = rng::randomReal();
        auto it = std::lower_bound(pop.begin(), pop.end(), threshold,
        [](const Candidate& sol, double threshold)
        {
            return sol.selection_cdf < threshold;
        });

        return (it != pop.end()) ? *it : pop.back();
    }

    template<typename geneType>
    inline typename GA<geneType>::Candidate GA<geneType>::sogaTournamentSelect(const Population& pop, size_t tourney_size)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
        assert(tourney_size > 1);

        /* Randomly pick tourney_size candidates. Indices may repeat. */
        std::vector<size_t> indices;
        indices.reserve(tourney_size);
        for (size_t i = 0; i < tourney_size; i++)
        {
            indices.push_back(rng::randomIdx(pop.size()));
        }

        /* Find the best of the picked candidates. */
        size_t idx = *std::max_element(indices.begin(), indices.end(),
        [&pop](size_t lidx, size_t ridx)
        {
            return pop[lidx].fitness[0] < pop[ridx].fitness[0];
        });

        return pop[idx];
    }

    template<typename geneType>
    inline typename GA<geneType>::Candidate GA<geneType>::sogaSelect(const Population& pop) const
    {
        switch (selection_method_)
        {
            case SogaSelection::tournament:
                return sogaTournamentSelect(pop, tournament_size_);
            case SogaSelection::roulette:
                [[fallthrough]];
            case SogaSelection::rank:
                [[fallthrough]];
            case SogaSelection::sigma:
                [[fallthrough]];
            case SogaSelection::boltzmann:
                return sogaWeightProportionalSelect(pop);
            case SogaSelection::custom:
                return customSelection(pop);
            default:
                assert(false);	/* Invalid selection method. Shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline typename GA<geneType>::Population GA<geneType>::updateSogaPopulation(Population& old_pop, CandidateVec& children) const
    {
        assert(old_pop.size() == population_size_);
        assert(!children.empty());
        assert(std::all_of(old_pop.begin(), old_pop.end(), [](const Candidate& sol) { return sol.is_evaluated; }));
        assert(std::all_of(children.begin(), children.end(), [](const Candidate& sol) { return sol.is_evaluated; }));

        old_pop.insert(old_pop.end(), make_move_iterator(children.begin()), make_move_iterator(children.end()));
        std::partial_sort(old_pop.begin(), old_pop.begin() + population_size_, old_pop.end(),
        [](const Candidate& lhs, const Candidate& rhs)
        {
            return lhs.fitness[0] > rhs.fitness[0];
        });
        old_pop.resize(population_size_);

        return old_pop;
    }
    
    template<typename geneType>
    inline std::vector<std::vector<size_t>> GA<geneType>::nonDominatedSort(Population& pop)
    {
        using namespace std;

        /* Calc the number of candidates which dominate each candidate, and the indices of the candidates it dominates. */
        vector<size_t> dom_count(pop.size(), 0);
        vector<vector<size_t>> dom_list(pop.size());

        for (size_t i = 0; i < pop.size(); i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                if (detail::paretoCompare(pop[j].fitness, pop[i].fitness))
                {
                    dom_count[j]++;
                    dom_list[i].push_back(j);
                }
                else if (detail::paretoCompare(pop[i].fitness, pop[j].fitness))
                {
                    dom_count[i]++;
                    dom_list[j].push_back(i);
                }
            }
        }

        /* Find the indices of all non-dominated candidates (first/best pareto front). */
        vector<size_t> front;
        for (size_t i = 0; i < pop.size(); i++)
        {
            if (dom_count[i] == 0)
            {
                front.push_back(i);
                pop[i].rank = 0;
            }
        }
        /* Find all the other pareto fronts. */
        vector<vector<size_t>> pareto_fronts;
        size_t front_idx = 1;
        while (!front.empty())
        {
            /* "Remove" the current front and find the next one. */
            vector<size_t> next_front;
            for (const auto& i : front)
            {
                for (const auto& j : dom_list[i])
                {
                    /* j belongs to the next front if it's domination count will become 0. */
                    if (--dom_count[j] == 0)
                    {
                        next_front.push_back(j);
                        pop[j].rank = front_idx;
                    }
                }
            }
            pareto_fronts.push_back(front);
            front = next_front;
            front_idx++;
        }

        return pareto_fronts;
    }

    template<typename geneType>
    inline void GA<geneType>::calcCrowdingDistances(Population& pop, std::vector<std::vector<size_t>>& pfronts)
    {
        using namespace std;
        assert(!pop.empty());

        for (const auto& pfront : pfronts)
        {
            for (const auto& idx : pfront)
            {
                pop[idx].distance = 0.0;
            }
        }

        for_each(execution::par_unseq, pfronts.begin(), pfronts.end(),
        [&pop](vector<size_t>& pfront)
        {
            /* Calc the distances in each fitness dimension. */
            for (size_t d = 0; d < pop[0].fitness.size(); d++)
            {
                sort(pfront.begin(), pfront.end(),
                [&pop, &d](size_t lidx, size_t ridx)
                {
                    return pop[lidx].fitness[d] < pop[ridx].fitness[d];
                });

                /* Calc the crowding distance for each solution. */
                double finterval = pop[pfront.back()].fitness[d] - pop[pfront.front()].fitness[d];
                finterval = max(finterval, 1E-6);

                pop[pfront.front()].distance = numeric_limits<double>::infinity();
                pop[pfront.back()].distance = numeric_limits<double>::infinity();
                for (size_t i = 1; i < pfront.size() - 1; i++)
                {
                    pop[pfront[i]].distance += (pop[pfront[i + 1]].fitness[d] - pop[pfront[i - 1]].fitness[d]) / finterval;
                }
            }
        });
    }

    template<typename geneType>
    inline bool GA<geneType>::crowdedCompare(const Candidate& lhs, const Candidate& rhs)
    {
        if (rhs.rank > lhs.rank) return true;
        else if (lhs.rank == rhs.rank) return lhs.distance > rhs.distance;
        else return false;
    }

    template<typename geneType>
    inline typename GA<geneType>::Candidate GA<geneType>::nsga2Select(const Population& pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop.size());
        size_t idx2 = rng::randomIdx(pop.size());

        return crowdedCompare(pop[idx1], pop[idx2]) ? pop[idx1] : pop[idx2];
    }

    template<typename geneType>
    inline typename GA<geneType>::Population GA<geneType>::updateNsga2Population(Population& old_pop, CandidateVec& children) const
    {
        using namespace std;
        assert(old_pop.size() == population_size_);
        assert(!children.empty());
        assert(all_of(old_pop.begin(), old_pop.end(), [](const Candidate& sol) { return sol.is_evaluated; }));
        assert(all_of(children.begin(), children.end(), [](const Candidate& sol) { return sol.is_evaluated; }));

        Population new_pop;
        new_pop.reserve(population_size_);

        old_pop.insert(old_pop.end(), make_move_iterator(children.begin()), make_move_iterator(children.end()));
        vector<vector<size_t>> pareto_fronts = nonDominatedSort(old_pop);
        calcCrowdingDistances(old_pop, pareto_fronts);

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        while (new_pop.size() + pareto_fronts[front_idx].size() <= population_size_)
        {
            for (const auto& idx : pareto_fronts[front_idx])
            {
                new_pop.push_back(move(old_pop[idx]));
            }
            front_idx++;
        }

        /* Add the remaining candidates from the partial front if there is one. */
        if (new_pop.size() != population_size_)
        {
            vector<size_t> added_indices(population_size_ - new_pop.size());	/* For updating the crowding distances in this front. */
            iota(added_indices.begin(), added_indices.end(), new_pop.size());

            vector<size_t> partial_front = pareto_fronts[front_idx];

            sort(partial_front.begin(), partial_front.end(),
            [&old_pop](size_t lidx, size_t ridx)
            {
                return crowdedCompare(old_pop[lidx], old_pop[ridx]);
            });

            for (const auto& idx : partial_front)
            {
                new_pop.push_back(move(old_pop[idx]));
                if (new_pop.size() == population_size_) break;
            }

            vector<vector<size_t>> temp = { added_indices };
            calcCrowdingDistances(new_pop, temp);
        }

        return new_pop;
    }

    template<typename geneType>
    inline void GA<geneType>::updateIdealPoint(const Population& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [this](const Candidate& sol) { return sol.fitness.size() == ideal_point_.size(); }));

        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < ideal_point_.size(); i++)
            {
                ideal_point_[i] = std::max(ideal_point_[i], sol.fitness[i]);
            }
        }
    }

    template<typename geneType>
    inline void GA<geneType>::updateNadirPoint(const Population& pop)
    {
        using namespace std;
        assert(!pop.empty());
        assert(all_of(pop.begin(), pop.end(), [this, &pop](const Candidate& sol) { return sol.fitness.size() == nadir_point_.size(); }));

        /* Identify/update extreme points for each objective axis. */
        for (size_t i = 0; i < nadir_point_.size(); i++)
        {
            vector<double> weights(nadir_point_.size(), 1E-6);
            weights[i] = 1.0;

            /* Find the solution or extreme point with the lowest Chebysev distance to the objective axis. */
            double dmin = numeric_limits<double>::max();
            vector<double> argmin;
            for (const auto& sol : pop)
            {
                double d = detail::ASF(sol.fitness, ideal_point_, weights);

                if (d < dmin)
                {
                    dmin = d;
                    argmin = sol.fitness;
                }
            }

            /* There are no extreme points yet in the first generation. */
            if (generation_cntr_ != 0)
            {
                for (const auto& extreme_point : extreme_points_)
                {
                    double d = detail::ASF(extreme_point, ideal_point_, weights);
                    if (d < dmin)
                    {
                        dmin = d;
                        argmin = extreme_point;
                    }
                }
            }

            extreme_points_[i] = argmin;
        }

        /* Find minimum of extreme points along each objective (nadir point). */
        for (size_t i = 0; i < nadir_point_.size(); i++)
        {
            nadir_point_[i] = extreme_points_[0][i];
            for (size_t j = 1; j < extreme_points_.size(); j++)
            {
                nadir_point_[i] = min(nadir_point_[i], extreme_points_[j][i]);
            }
        }
    }

    template<typename geneType>
    inline void GA<geneType>::associatePopToRefs(Population& pop, const std::vector<std::vector<double>>& ref_points)
    {
        using namespace std;
        assert(!pop.empty());
        assert(all_of(pop.begin(), pop.end(), [&pop](const Candidate& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        updateIdealPoint(pop);
        updateNadirPoint(pop);

        vector<vector<double>> fnorms(pop.size(), vector<double>(pop[0].fitness.size(), 0.0));	/* Don't change the actual fitness values. */

        transform(execution::par_unseq, pop.begin(), pop.end(), fnorms.begin(), fnorms.begin(),
        [this](const Candidate& sol, vector<double>& fnorm) -> vector<double>
        {
            for (size_t i = 0; i < sol.fitness.size(); i++)
            {
                fnorm[i] = sol.fitness[i] - ideal_point_[i];
                fnorm[i] /= min(nadir_point_[i] - ideal_point_[i], -1E-6);
            }

            return fnorm;
        });

        /* Associate each candidate with the closest reference point. */
        transform(execution::par_unseq, pop.begin(), pop.end(), fnorms.begin(), pop.begin(),
        [&ref_points](Candidate& sol, const vector<double>& f) -> Candidate
        {
            tie(sol.ref_idx, sol.distance) = detail::findClosestRef(ref_points, f);
            return sol;
        });
    }

    template<typename geneType>
    inline std::vector<size_t> GA<geneType>::calcNicheCounts(Population& pop, const std::vector<std::vector<double>>& ref_points)
    {
        std::vector<size_t> niche_counts(ref_points.size(), 0U);
        for (const auto& sol : pop)
        {
            niche_counts[sol.ref_idx]++;
        }

        /* Assign the niche counts to the candidates too. */
        for (auto& sol : pop)
        {
            sol.niche_count = niche_counts[sol.ref_idx];
        }

        return niche_counts;
    }

    template<typename geneType>
    inline bool GA<geneType>::nichedCompare(const Candidate& lhs, const Candidate& rhs)
    {
        if (rhs.rank > lhs.rank) return true;
        else if (lhs.rank == rhs.rank) return lhs.niche_count < rhs.niche_count;
        else if (lhs.rank == rhs.rank && lhs.niche_count == rhs.niche_count) return lhs.distance < rhs.distance;
        else return false;
    }

    template<typename geneType>
    inline typename GA<geneType>::Candidate GA<geneType>::nsga3Select(const Population& pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop.size());
        size_t idx2 = rng::randomIdx(pop.size());

        return nichedCompare(pop[idx1], pop[idx2]) ? pop[idx1] : pop[idx2];
    }

    template<typename geneType>
    inline typename GA<geneType>::Population GA<geneType>::updateNsga3Population(Population& old_pop, CandidateVec& children)
    {
        using namespace std;
        assert(old_pop.size() == population_size_);
        assert(!children.empty());
        assert(all_of(old_pop.begin(), old_pop.end(), [](const Candidate& sol) { return sol.is_evaluated; }));
        assert(all_of(children.begin(), children.end(), [](const Candidate& sol) { return sol.is_evaluated; }));

        Population new_pop;
        new_pop.reserve(population_size_);

        old_pop.insert(old_pop.end(), make_move_iterator(children.begin()), make_move_iterator(children.end()));
        vector<vector<size_t>> pareto_fronts = nonDominatedSort(old_pop);
        associatePopToRefs(old_pop, ref_points_);

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        while (new_pop.size() + pareto_fronts[front_idx].size() <= population_size_)
        {
            for (const auto& idx : pareto_fronts[front_idx])
            {
                new_pop.push_back(move(old_pop[idx]));
            }
            front_idx++;
        }
        vector<size_t> niche_counts = calcNicheCounts(new_pop, ref_points_);

        /* Add remaining candidates from the partial front if there is one. */
        vector<size_t> partial_front = pareto_fronts[front_idx];
        while (new_pop.size() != population_size_)
        {
            /* Find the lowest niche count in the partial front. */
            size_t min_count = population_size_;
            for (const auto& idx : partial_front)
            {
                min_count = min(min_count, niche_counts[old_pop[idx].ref_idx]);
            }

            /* Find the reference points with minimal niche counts, and pick one. */
            vector<size_t> refs = {};
            for (const auto& idx : partial_front)
            {
                size_t ref = old_pop[idx].ref_idx;
                if (niche_counts[ref] == min_count && find(refs.begin(), refs.end(), ref) == refs.end())
                {
                    refs.push_back(ref);
                }
            }
            size_t ref = refs[rng::randomIdx(refs.size())];

            /* Find the idx of the closest sol in the partial front associated with this ref point. */
            size_t sol_idx = partial_front[0];
            double min_distance = numeric_limits<double>::infinity();
            for (const auto& idx : partial_front)
            {
                if (old_pop[idx].ref_idx == ref && old_pop[idx].distance < min_distance)
                {
                    min_distance = old_pop[idx].distance;
                    sol_idx = idx;
                }
            }

            /* Move this candidate to new_pop and increment the associated niche count. */
            new_pop.push_back(move(old_pop[sol_idx]));
            partial_front.erase(remove(partial_front.begin(), partial_front.end(), sol_idx), partial_front.end());

            niche_counts[ref]++;
            for (auto& sol : new_pop)
            {
                if (sol.ref_idx == ref) sol.niche_count++;
            }
        }

        return new_pop;
    }

    template<typename geneType>
    inline typename GA<geneType>::CandidateVec GA<geneType>::findParetoFront1D(const Population& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1; }));

        CandidateVec optimal_sols;

        double fmax = fitnessMax(pop)[0]; /* There might be multiple solutions with this max fitness value. */
        for (const auto& sol : pop)
        {
            if (sol.fitness[0] == fmax) optimal_sols.push_back(sol);
        }

        return optimal_sols;
    }

    template<typename geneType>
    inline typename GA<geneType>::CandidateVec GA<geneType>::findParetoFrontKung(const Population& pop)
    {
        /* See: Kung et al. "On finding the maxima of a set of vectors." Journal of the ACM (JACM) 22.4 (1975): 469-476.*/
        /* Doesn't work for d = 1 (single-objective optimization). */

        using namespace std;
        using iter = vector<size_t>::iterator;

        assert(!pop.empty());
        assert(all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return !sol.fitness.empty(); }));
        assert(all_of(pop.begin(), pop.end(), [&pop](const Candidate& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        size_t dim = pop[0].fitness.size();	/* The number of objectives. */

        /* Find the indices of pareto optimal solutions in the population (Kung's algorithm) assuming fitness maximization. */
        function<vector<size_t>(iter, iter)> pfront = [&pfront, &pop, &dim](iter first, iter last) -> vector<size_t>
        {
            if (distance(first, last) == 1) return { *first };

            vector<size_t> R = pfront(first, first + distance(first, last) / 2);	/* Top half. */
            vector<size_t> S = pfront(first + distance(first, last) / 2, last);		/* Bottom half. */

            /* T = Find all non-dominated elements of the bottom half. */
            vector<size_t> T;
            for (const auto& s : S)
            {
                /* Check if s is dominated by any solution in R. */
                bool is_dominated = false;
                for (const auto& r : R)
                {
                    /* Pareto compare s and r. */
                    /* The first dimension (d = 0) of the fitness vectors doesn't need to be compared since the pop is already sorted. */
                    for (size_t d = 1; d < dim; d++)
                    {
                        if (pop[s].fitness[d] > pop[r].fitness[d])
                        {
                            is_dominated = false;
                            break;
                        }
                        if (pop[s].fitness[d] < pop[r].fitness[d]) is_dominated = true;
                    }
                    if (is_dominated) break;
                }
                if (!is_dominated) T.push_back(s);
            }
            R.insert(R.end(), T.begin(), T.end());

            return R;
        };

        /* Find the indices of the pareto optimal candidates. */
        vector<size_t> indices(pop.size());
        iota(indices.begin(), indices.end(), 0U);

        /* Sort the pop indices into descending order based on first fitness value (needed for Kung's). */
        sort(indices.begin(), indices.end(), [&pop](size_t lidx, size_t ridx) { return pop[lidx].fitness[0] > pop[ridx].fitness[0]; });
        indices = pfront(indices.begin(), indices.end());

        CandidateVec optimal_sols;
        optimal_sols.reserve(indices.size());
        for (const auto& idx : indices)
        {
            optimal_sols.push_back(pop[idx]);
        }

        return optimal_sols;
    }

    template<typename geneType>
    inline std::vector<double> GA<geneType>::fitnessMin(const Population& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return !sol.fitness.empty(); }));
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const Candidate& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        std::vector<double> fmin = pop[0].fitness;

        for (size_t i = 1; i < pop.size(); i++)
        {
            for (size_t j = 0; j < fmin.size(); j++)
            {
                fmin[j] = std::min(fmin[j], pop[i].fitness[j]);
            }
        }

        return fmin;
    }

    template<typename geneType>
    inline std::vector<double> GA<geneType>::fitnessMax(const Population& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return !sol.fitness.empty(); }));
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const Candidate& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        std::vector<double> fmax = pop[0].fitness;

        for (size_t i = 1; i < pop.size(); i++)
        {
            for (size_t j = 0; j < fmax.size(); j++)
            {
                fmax[j] = std::max(fmax[j], pop[i].fitness[j]);
            }
        }

        return fmax;
    }

    template<typename geneType>
    inline double GA<geneType>::fitnessMean(const Population& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return !sol.fitness.empty(); }));

        return std::accumulate(pop.begin(), pop.end(), 0.0,
        [&pop](double sum, const Candidate& sol)
        {
            return sum + sol.fitness[0] / pop.size();
        });
    }

    template<typename geneType>
    inline double GA<geneType>::fitnessSD(const Population& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return !sol.fitness.empty(); }));

        if (pop.size() == 1) return 0.0;

        double mean = fitnessMean(pop);
        long double variance = std::accumulate(pop.begin(), pop.end(), 0.0L,
        [&pop, mean](long double sum, const Candidate& sol)
        {
            return sum + std::pow(sol.fitness[0] - mean, 2) / (pop.size() - 1.0);
        });

        return double(std::sqrt(variance));
    }

} // namespace genetic_algorithm

#endif // !GA_BASE_GA_H