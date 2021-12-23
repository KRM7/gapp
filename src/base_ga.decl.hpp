#ifndef GA_BASE_GA_DECL_HPP
#define GA_BASE_GA_DECL_HPP


#include <algorithm>
#include <vector>
#include <unordered_set>
#include <utility>
#include <functional>
#include <atomic>
#include <cstddef>

#include <memory>
#include "candidate.hpp"


namespace genetic_algorithm
{
    namespace crossover
    {
        template<regular_hashable GeneType>
        class Crossover;
    }
    namespace mutation
    {
        template<regular_hashable GeneType>
        class Mutation;
    }

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
            std::vector<double> fitness_mean;   /**< The mean fitness values of each generation. */
            std::vector<double> fitness_sd;     /**< The standard deviation of the fitness values of each generation. */
            std::vector<double> fitness_min;    /**< The lowest fitness value in each generation. */
            std::vector<double> fitness_max;    /**< The highest fitness value in each generation. */

            void clear() noexcept;
            void reserve(size_t new_capacity);
            void add(double mean, double sd, double min, double max);
        };

        using GeneType = geneType;
        using Candidate = Candidate<GeneType>;  /**< The candidates used in the algorithm, each representing a solution to the problem. */

        using Chromosome = std::vector<GeneType>;                                           /**< . */
        using CandidatePair = CandidatePair<GeneType>;                                      /**< . */
        using CandidateVec = std::vector<Candidate>;                                        /**< . */
        using CandidateSet = std::unordered_set<Candidate, CandidateHasher<GeneType>>;      /**< . */
        using Population = std::vector<Candidate>;                                          /**< . */

        using fitnessFunction_t = std::function<std::vector<double>(const Chromosome&)>;    /**< The type of the fitness function. */
        using selectionFunction_t = std::function<Candidate(const Population&)>;            /**< The type of the selection function. */
        using mutationFunction_t = std::function<void(Candidate&, double)>;                 /**< The type of the mutation function. */
        using repairFunction_t = std::function<Chromosome(const Chromosome&)>;              /**< The type of the repair function. */
        using callbackFunction_t = std::function<void(const GA*)>;

        /**
        * The type of the genetic algorithm used, depending on the problem type (single- or multi-objective optimization). \n
        * Set the mode used with @ref mode.
        */
        enum class Mode
        {
            single_objective,           /**< Simple single-objective genetic algorithm. */
            multi_objective_sorting,    /**< Non-dominated sorting genetic algorithm (NSGA-II) for multi-objective optimization. */
            multi_objective_decomp      /**< NSGA-III algorithm for many-objective optimization. */
        };

        /**
        * The possible stop conditions used in the algorithm. The algorithm always stops when @ref max_gen has been reached,
        * regardless of the stop condition selected. \n
        * Some of the stop condition do not work for multi-objective problems (fitness_mean_stall and fitness_best_stall).
        * Choose the stop condition with @ref stop_condition. \n
        */
        enum class StopCondition
        {
            max_gen,               /**< Only stop when @ref max_gen is reached. */
            fitness_value,         /**< Stop when a solution was found which dominates a reference fitness value. @see fitness_threshold */
            fitness_evals,         /**< Stop when the fitness function has been evaluated a set number of times. @see max_fitness_evals */
            fitness_mean_stall,    /**< Stop when the mean fitness of the population doesn't improve at least @ref stall_threshold over @ref stall_gen_count. */
            fitness_best_stall     /**< Stop when the highest fitness of the population doesn't improve at least @ref stall_threshold over @ref stall_gen_count. */
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
            roulette,          /**< Standard roulette selection adapted to also work with negative fitness values. No parameters. */
            rank,              /**< Standard rank selection. @see rank_sel_min_w_, @see rank_sel_max_w_ */
            tournament,        /**< Standard tournament selection. @see tournament_size */
            sigma,             /**< Sigma fitness scaling. @see sigma_scale */
            boltzmann,         /**< Standard Boltzmann selection. @see boltzmann_temps */
            custom             /**< A user defined function is used to compute the selection probabilities. @see customCalcWeights */
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
        [[maybe_unused]] const CandidateVec& run();


        /** @returns A vector of the pareto optimal solutions found while running the algorithm. */
        [[nodiscard]] const CandidateVec& solutions() const;

        /** @returns The number of fitness evaluations performed while running the algorithm. */
        [[nodiscard]] size_t num_fitness_evals() const;

        /** @returns The current value of the generation counter. */
        [[nodiscard]] size_t generation_cntr() const;

        /** @returns The population of the final generation in the algorithm. */
        [[nodiscard]] const Population& population() const;

        /** @returns A History object containing stats from each generation of the single objective genetic algorithm. */
        [[nodiscard]] const History& soga_history() const;

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

        /* CROSSOVER */
        template<typename CrossoverType>
        //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>> && std::copy_constructible<CrossoverType>
        void crossover_method(const CrossoverType& f);

        void crossover_method(std::unique_ptr<crossover::Crossover<GeneType>>&& f);

        template<typename CrossoverType = crossover::Crossover<GeneType>>
        //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>>
        CrossoverType& crossover_method() const;


        /* MUTATION */
        template<typename MutationType>
        //requires std::derived_from<MutationType, mutation::Mutation<GeneType>> && std::copy_constructible<MutationType>
        void mutation_method(const MutationType& f);

        void mutation_method(std::unique_ptr<mutation::Mutation<GeneType>>&& f);

        template<typename MutationType = mutation::Mutation<GeneType>>
        //requires std::derived_from<MutationType, mutation::Crossover<GeneType>>
        MutationType& mutation_method() const;

    protected:

        Population population_;
        size_t generation_cntr_ = 0;
        size_t num_objectives_ = 0;        /* Determined from the fitness function. */

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

        std::unique_ptr<crossover::Crossover<GeneType>> crossover_;
        std::unique_ptr<mutation::Mutation<geneType>> mutation_;


        /* General functions for the genetic algorithms. */

        void init();
        virtual Candidate generateCandidate() const = 0;
        Population generateInitialPopulation() const;
        void evaluate(Population& pop);
        void updateOptimalSolutions(CandidateVec& optimal_sols, const Population& pop) const;
        void prepSelections(Population& pop) const;
        Candidate select(const Population& pop) const;
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

#endif // !GA_BASE_GA_DECL_HPP