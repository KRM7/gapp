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
* Implementations of genetic algorithms with binary, real, permutational, and integer encodings.
* 
* @file genetic_algorithm.h
*/

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <random>
#include <functional>
#include <numeric>
#include <limits>
#include <execution>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cassert>
#include <stdexcept>

/** Genetic algorithms and random number generation. */
namespace genetic_algorithm
{
    /** Contains the PRNG classes and functions for generating random numbers. */
    namespace rng
    {
        /**
        * Splitmix64 PRNG adapted from https://prng.di.unimi.it/splitmix64.c \n
        * Only used for the initialization of the other PRNGs.
        */
        class splitmix64
        {
        public:
            using result_type = uint_fast64_t;
            using state_type = uint_fast64_t;

            splitmix64(state_type seed) : state(seed) {}

            /* Generate the next random number. */
            result_type operator()() noexcept
            {
                state += 0x9e3779b97f4a7c15;
                result_type z = state;
                z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
                z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
                return z ^ (z >> 31);
            }

        private:
            state_type state;
        };

        /**
        * xoroshiro128+ PRNG adapted from https://prng.di.unimi.it/xoroshiro128plus.c \n
        * Works with the standard library distributions, qualifies as std::uniform_random_bit_generator
        */
        class xoroshiro128p
        {
        public:
            using result_type = uint_fast64_t;
            using state_type = uint_fast64_t;

            xoroshiro128p(uint_fast64_t seed)
            {
                splitmix64 seed_seq_gen(seed);
                state[0] = seed_seq_gen();
                state[1] = seed_seq_gen();
            }

            result_type operator()() noexcept
            {
                state_type s0 = state[0];
                state_type s1 = state[1];
                result_type result = s0 + s1;

                s1 ^= s0;
                state[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16);
                state[1] = rotl(s1, 37);

                return result;
            }

            static constexpr result_type min() noexcept
            {
                return std::numeric_limits<result_type>::lowest();
            }
            static constexpr result_type max() noexcept
            {
                return std::numeric_limits<result_type>::max();
            }

        private:
            state_type state[2];

            static state_type rotl(state_type x, int k) noexcept
            {
                return (x << k) | (x >> (64 - k));
            }
        };

        /** The PRNG used in the genetic algorithm. */
        using PRNG = xoroshiro128p;

        /**
        * Generates a random double on the interval [l_bound, u_bound).
        * 
        * @param l_bound The lower bound of the interval.
        * @param u_bound The upper bound of the interval.
        * @return The generated random number.
        */
        inline double generateRandomDouble(double l_bound = 0.0, double u_bound = 1.0)
        {
            assert(l_bound <= u_bound);

            static thread_local PRNG engine{ std::random_device{}() };
            std::uniform_real_distribution<double> distribution{ l_bound, u_bound };

            return distribution(engine);
        }

        /**
        * Generates a random integer on the closed interval [l_bound, u_bound].
        *
        * @tparam T The generated integer type.
        * @param l_bound The lower bound of the interval.
        * @param u_bound The upper bound of the interval.
        * @return The generated random number.
        */
        template<typename T>
        inline T generateRandomInt(T l_bound, T u_bound)
        {
            assert(l_bound <= u_bound);

            static thread_local PRNG engine{ std::random_device{}() };
            std::uniform_int_distribution<T> distribution{ l_bound, u_bound };

            return distribution(engine);
        }

        /**
        * Generates a random unsigned integer on the closed interval [0, c_size-1]. \n
        * Used to generate a random index for containers, with c_size being the size of the container.
        *
        * @param c_size The size of the container.
        * @return The generated random number.
        */
        inline size_t generateRandomIdx(size_t c_size)
        {
            return generateRandomInt(size_t{ 0 }, c_size - 1);
        }

        /**
        * Generates a random boolean value.
        *
        * @return The generated random boolean.
        */
        inline bool generateRandomBool()
        {
            return bool(generateRandomInt(size_t{ 0 }, size_t{ 1 }));
        }

        /**
        * Generates a random double from a normal distribution.
        *
        * @param mean The mean of the normal distribution.
        * @param SD The standard deviation of the normal distribution.
        * @return The generated random number.
        */
        inline double generateRandomNorm(double mean = 0.0, double SD = 1.0)
        {
            assert(SD > 0.0);

            static thread_local PRNG engine{ std::random_device{}() };
            std::normal_distribution<double> distribution{ mean, SD };

            return distribution(engine);
        }

        /**
        * Sample a point from a uniform distribution on a unit simplex in @p dim dimensions. \n
        * 
        * @param dim The number of elements in the generated vector.
        * @returns The generated random point/vector.
        */
        inline std::vector<double> generateRandomSimplexPoint(size_t dim)
        {
            assert(dim > 0);

            static thread_local std::minstd_rand0 engine{ std::random_device{}() };
            std::uniform_real_distribution<double> distribution{ 0.0, 1.0 };

            std::vector<double> vec;
            vec.reserve(dim);

            double sum = 0.0;
            for (size_t i = 0; i < dim; i++)
            {
                vec.push_back(-std::log(distribution(engine)));
                sum += vec.back();
            }
            for (size_t i = 0; i < dim; i++)
            {
                vec[i] /= sum;
            }

            return vec;
        }

    } // namespace rng

    namespace detail
    {
        /* Return true if lhs is dominated by rhs (lhs < rhs) assuming maximization. */
        inline bool paretoCompare(const std::vector<double>& lhs, const std::vector<double>& rhs)
        {
            assert(lhs.size() == rhs.size());

            bool has_lower = false;
            for (size_t i = 0; i < lhs.size(); i++)
            {
                if (lhs[i] > rhs[i]) return false;
                if (lhs[i] < rhs[i]) has_lower = true;
            }
            return has_lower;
        }

        /* Calculate the square of the Euclidean distance between the vectors v1 and v2. */
        inline double euclideanDistanceSq(const std::vector<double>& v1, const std::vector<double>& v2)
        {
            assert(v1.size() == v2.size());

            double d = 0.0;
            for (size_t i = 0; i < v1.size(); i++)
            {
                d += std::pow(v1[i] - v2[i], 2);
            }

            return d;
        }

        /* Generate n reference points in dim dimensions for the NSGA-III algorithm. */
        inline std::vector<std::vector<double>> generateRefPoints(size_t n, size_t dim)
        {
            using namespace std;
            assert(n > 0);
            assert(dim > 1);

            /* Generate random reference point candidates. */
            size_t k = max(size_t{ 10 }, 2 * dim);
            vector<vector<double>> candidates(k * n - 1);
            generate(candidates.begin(), candidates.end(), [&dim]() { return rng::generateRandomSimplexPoint(dim); });

            vector<vector<double>> refs;
            refs.reserve(n);
            /* The first ref point is random. */
            refs.push_back(rng::generateRandomSimplexPoint(dim));

            vector<double> min_distances(candidates.size(), numeric_limits<double>::infinity());
            while (refs.size() < n)
            {
                /* Calc the distance of each candidate to the closest ref point. */
                transform(execution::par_unseq, candidates.begin(), candidates.end(), min_distances.begin(), min_distances.begin(),
                [&refs](const vector<double>& candidate, double dmin)
                {
                    double d = euclideanDistanceSq(candidate, refs.back());
                    return min(dmin, d);
                });

                /* Add the candidate with highest min_distance to the refs. */
                size_t idx = size_t(max_element(min_distances.begin(), min_distances.end()) - min_distances.begin());
                refs.push_back(move(candidates[idx]));

                /* Remove the added candidate and the corresponding min_distance. */
                swap(candidates[idx], candidates.back());
                candidates.pop_back();
                swap(min_distances[idx], min_distances.back());
                min_distances.pop_back();
            }

            return refs;
        }

        /* Calculate the square of the perpendicular distance between the line ref and the point p. */
        inline double perpendicularDistanceSq(const std::vector<double>& ref, const std::vector<double>& p)
        {
            assert(ref.size() == p.size());

            double num = 0.0, den = 0.0;
            for (size_t i = 0; i < ref.size(); i++)
            {
                num += ref[i] * p[i];
                den += ref[i] * ref[i];
            }
            double k = num / den;

            double dist = 0.0;
            for (size_t i = 0; i < ref.size(); i++)
            {
                dist += std::pow(p[i] - k * ref[i], 2);
            }

            return dist;
        }

        /* Find the index and distance of the closest reference line to the point p. */
        inline std::pair<size_t, double> findClosestRef(const std::vector<std::vector<double>>& refs, const std::vector<double>& p)
        {
            size_t argmin = 0;
            double dmin = perpendicularDistanceSq(refs[0], p);
            for (size_t i = 1; i < refs.size(); i++)
            {
                double d = perpendicularDistanceSq(refs[i], p);
                if (d < dmin)
                {
                    dmin = d;
                    argmin = i;
                }
            }

            return std::make_pair(argmin, dmin);
        }

        /* Achievement scalarization function. */
        inline double ASF(const std::vector<double>& f, const std::vector<double>& z, const std::vector<double>& w)
        {
            assert(!f.empty());
            assert(f.size() == z.size() && f.size() == w.size());

            double dmax = std::abs(f[0] - z[0]) / w[0];
            for (size_t j = 1; j < f.size(); j++)
            {
                dmax = std::max(dmax, std::abs(f[j] - z[j]) / w[j]);
            }

            return dmax;
        }

    } // namespace detail

    /**
    * Abstract base GA class.
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

            void clear() noexcept
            {
                fitness_mean.clear();
                fitness_sd.clear();
                fitness_min.clear();
                fitness_max.clear();
            }

            void reserve(size_t new_capacity)
            {
                fitness_mean.reserve(new_capacity);
                fitness_sd.reserve(new_capacity);
                fitness_min.reserve(new_capacity);
                fitness_max.reserve(new_capacity);
            }

            void add(double mean, double sd, double min, double max)
            {
                fitness_mean.push_back(mean);
                fitness_sd.push_back(sd);
                fitness_min.push_back(min);
                fitness_max.push_back(max);
            }
        };

        /** The candidates used in the algorithm, each representing a solution to the problem. */
        struct Candidate
        {
            std::vector<geneType> chromosome;	/**< The chromosome encoding the solution. */
            std::vector<double> fitness;		/**< The fitness values of the candidate solution. */

            double selection_pdf = 0.0;			/**< The probability of selecting the candidate (SOGA). */
            double selection_cdf = 0.0;			/**< The value of the cumulative distribution function for the candidate (SOGA). */

            size_t rank = 0;					/**< Non-domination rank. (used in both the NSGA-II and NSGA-III) */
            double distance = 0.0;				/**< Crowding distance (NSGA-II), or the distance to closest reference point (NSGA-III). */
            size_t ref_idx = 0;					/**< Index of the associated reference point (NSGA-III). */
            size_t niche_count = 0;				/**< Number of candidates associated with the same reference point as this candidate (NSGA-III). */

            bool is_evaluated = false;			/**< False if the candidate's fitness value needs to be computed. */

            Candidate();
            Candidate(const std::vector<geneType>& chrom) : chromosome(chrom) {}
            Candidate(std::vector<geneType>&& chrom) noexcept : chromosome(std::move(chrom)) {}

            bool operator==(const Candidate& rhs) const
            {
                return (this->chromosome == rhs.chromosome);
            }
            bool operator!=(const Candidate& rhs) const
            {
                return (this->chromosome != rhs.chromosome);
            }
        };

        /**
        * Hasher for the candidates so they can be stored in an unordered set. \n
        * A hash function must be defined for the geneType for this to work.
        */
        struct CandidateHasher
        {
            size_t operator()(const Candidate& c) const noexcept
            {
                size_t seed = c.chromosome.size();
                for (const auto& gene : c.chromosome)
                {
                    seed ^= std::hash<geneType>()(gene) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                }
                return seed;
            }
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
        * Used to eliminate pointless fitness evaluations.
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
        GA(size_t chrom_len, fitnessFunction_t fitness_function) : 
            chrom_len_(chrom_len),
            mutation_rate_(1.0 / chrom_len),
            fitnessFunction(fitness_function)
        {
            if (chrom_len == 0) throw std::invalid_argument("The chromosome length must be at least 1.");
            if (fitnessFunction == nullptr) throw std::invalid_argument("The fitness function is a nullptr.");
        };

        GA() = delete;
        GA(const GA&) = default;
        GA& operator=(const GA&) = default;
        GA(GA&&) = default;
        GA& operator=(GA&&) = default;
        virtual ~GA() = default;

        /**
        * Runs the genetic algorithm with the selected settings.
        * 
        * @returns The optimal solutions.
        */
        [[maybe_unused]] CandidateVec run()
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

        /**
        * @returns A vector of the pareto optimal solutions found while running the algorithm.
        */
        [[nodiscard]] CandidateVec solutions() const { return solutions_; }

        /**
        * @returns The number of fitness evaluations performed while running the algorithm.
        */
        [[nodiscard]] size_t num_fitness_evals() const { return static_cast<size_t>(num_fitness_evals_); }

        /**
        * @returns The current value of the generation counter.
        */
        [[nodiscard]] size_t generation_cntr() const { return generation_cntr_; }

        /**
        * @returns The population of the final generation in the algorithm.
        */
        [[nodiscard]] Population population() const { return population_; }

        /**
        * @returns A History object containing stats from each generation of the single objective genetic algorithm.
        */
        [[nodiscard]] History soga_history() const { return soga_history_; }


        /**
        * Set the type of the problem/genetic algorithm that will be used (single-/multi-objective).
        *
        * @param mode The algorithm type to use. @see Mode
        */
        void mode(Mode mode)
        {
            if (static_cast<size_t>(mode) > 2) throw std::invalid_argument("Invalid algorithm mode selected.");

            mode_ = mode;
        }
        [[nodiscard]] Mode mode() const { return mode_; }

        /**
        * Sets the length of the chromosomes (number of genes) of the Candidate solutions used in the algorithm to @p len. \n
        * The chromosome length must be at least 1.
        *
        * @param len The length of the chromosomes.
        */
        void chrom_len(size_t len)
        {
            if (len == 0) throw std::invalid_argument("The chromosome length must be at least 1.");

            chrom_len_ = len;
        }
        [[nodiscard]] size_t chrom_len() const { return chrom_len_; }

        /**
        * Sets the number of Candidates used in the population to @p size. \n
        * The population size must be at least 1.
        *
        * @param size The size of the populations.
        */
        void population_size(size_t size)
        {
            if (size == 0) throw std::invalid_argument("The population size must be at least 1.");

            population_size_ = size;
        }
        [[nodiscard]] size_t population_size() const { return population_size_; }

        /**
        * Sets the crossover rate used in the algorithm to @p pc. \n
        * The crossover rate must be on the closed interval [0.0, 1.0].
        *
        * @param pc The probability of crossover being performed on 2 selected individuals.
        */
        void crossover_rate(double pc)
        {
            if (!(0.0 <= pc && pc <= 1.0)) throw std::invalid_argument("The crossover probability must be in the range [0, 1].");

            crossover_rate_ = pc;
        }
        [[nodiscard]] double crossover_rate() const { return crossover_rate_; }

        /**
        * Sets the mutation rate used in the algorithm to @p pm. \n
        * The mutation rate must be on the closed interval [0.0, 1.0].
        *
        * @param pm The probability of mutating the children generated by the crossovers (per child or per gene based on encoding type).
        */
        void mutation_rate(double pm)
        {
            if (!(0.0 <= pm && pm <= 1.0)) throw std::invalid_argument("The mutation probability must be in the range [0, 1].");

            mutation_rate_ = pm;
        }
        [[nodiscard]] double mutation_rate() const { return mutation_rate_; }

        /**
        * Sets the selection method used in the single-objective algorithm to @p method. \n
        * The selection method set is ignored in the other algorithm types. @see Mode
        *
        * @param method The selection method used in the single-objective algorithm.
        */
        void selection_method(SogaSelection method)
        {
            if (static_cast<size_t>(method) > 5) throw std::invalid_argument("Invalid soga selection method selected.");

            selection_method_ = method;
        }
        [[nodiscard]] SogaSelection selection_method() const { return selection_method_; }
        /**
        * Sets the selection function used in the single-objective algorithm to @p f. \n
        * The selection function set is ignored in the other algorithm types. @see Mode
        *
        * @param f The selection function used in the single-objective algorithm.
        */
        void selection_method(selectionFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The selection function can't be a nullptr.");

            selection_method_ = SogaSelection::custom;
            customSelection = f;
        }

        /**
        * Sets the number of solutions picked for the tournaments to @p size if the tournament selection method is
        * selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The size of the tournaments must be at least 2.
        *
        * @param size The size of the tournaments during tournament selection.
        */
        void tournament_size(size_t size)
        {
            if (size < 2) throw std::invalid_argument("The tournament size must be at least 2.");

            tournament_size_ = size;
        }
        [[nodiscard]] size_t tournament_size() const { return tournament_size_; }

        /**
        * Sets the minimum and maximum selection weights used during rank selection if the rank selection method is
        * selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The @p min_weight must be on the closed interval [0.0, @p max_weight]. \n
        * The @p max_weight must be greater than @p min_weight.
        *
        * @param min_weight The selection weight assigned to the worst Candidate in the population.
        * @param max_weight The selection weight assigned to the best Candidate in the population.
        */
        void rank_sel_weights(double min_weight, double max_weight)
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
        [[nodiscard]] std::pair<double, double> rank_sel_weights() const { return { rank_sel_min_w_, rank_sel_max_w_ }; }

        /**
        * Sets the minimum and maximum temperature values used during boltzmann selection if the boltzmann
        * selection method is selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The minimum @p tmin temperature must be on the interval [0.1, @p tmax). \n
        * The maximum @p tmax temperature must be on the interval (@p tmin, DBL_MAX).
        *
        * @param tmin The minimum temperature (used in the last generation).
        * @param tmax The maximum temperature (used in the first generation).
        */
        void boltzmann_temps(double tmin, double tmax)
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
        [[nodiscard]] std::pair<double, double> boltzmann_temps() const { return { boltzmann_tmin_, boltzmann_tmax_ }; }

        /**
        * Sets the scaling parameter used during selection to @p scale if the sigma selection method is
        * selected for the single-objective algorithm. @see selection_method @see SogaSelection \n
        * The scaled fitness values are: f' = (f - f_mean) / (@p scale * f_sd) \n
        * The value of @p scale must be on the interval [1.0, DBL_MAX].
        *
        * @param scale The scaling parameter used during sigma selection.
        */
        void sigma_scale(double scale)
        {
            if (!(1.0 <= scale && scale <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Scale must be in the range [1.0, DBL_MAX].");
            }

            sigma_scale_ = scale;
        }
        [[nodiscard]] double sigma_scale() const { return sigma_scale_; }

        /**
        * Sets the stop condition used in the algorithm to @p condition. Some of the stop conditions
        * only work with the single-objective algorithm. @see StopCondition \n
        * The algorithm always stops when the set max_gen generation has been reached regardless of the stop condition
        * set. @see max_gen
        *
        * @param condition The stop condition used in the algorithm.
        */
        void stop_condition(StopCondition condition)
        {
            if (static_cast<size_t>(condition) > 4) throw std::invalid_argument("Invalid stop condition selected.");

            stop_condition_ = condition;
        }
        [[nodiscard]] StopCondition stop_condition() const { return stop_condition_; }

        /**
        * Sets the maximum number of generations the algorithm runs for to @p max_gen. The
        * algorithm will always stop when this generation has been reached regardless of what stop
        * condition was set, but it can stop earlier when another stop condition is selected.
        * @see stop_condition @see StopCondition \n
        * The value of @p max_gen must be at least 1.
        *
        * @param max_gen The maximum number of generations.
        */
        void max_gen(size_t max_gen)
        {
            if (max_gen == 0) throw std::invalid_argument("The maximum number of generations must be at least 1.");

            max_gen_ = max_gen;
        }
        [[nodiscard]] size_t max_gen() const { return max_gen_; }

        /**
        * Sets the maximum number of fitness evaluations the algorithm runs for to @p max_evals if
        * the fitness_evals stop condition is selected. @see stop_condition @see StopCondition \n
        * The algorithm may evaluate the fitness function more times than the maximum set since the
        * stop condition is only checked at the end of each generation. \n
        * The value of @p max_evals must be at least 1.
        *
        * @param max_evals
        */
        void max_fitness_evals(size_t max_evals)
        {
            if (max_evals == 0) throw std::invalid_argument("The maximum number of fitness evaluations must be at least 1.");

            max_fitness_evals_ = max_evals;
        }
        [[nodiscard]] size_t max_fitness_evals() const { return max_fitness_evals_; }

        /**
        * Sets the reference fitness value for the fitness_value stop condition to @p ref. \n
        * The algorithm will stop running if a solution has been found which dominates this reference point. \n
        * The size of the reference vector should be equal to the number of objectives. \n
        * @see stop_condition @see StopCondition
        *
        * @param ref The fitness reference used.
        */
        void fitness_threshold(std::vector<double> ref)
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
        [[nodiscard]] std::vector<double> fitness_threshold() const { return fitness_reference_; }

        /**
        * Sets the number of generations to look back when evaluating the stall stop conditions. \n
        * Only relevant for the single-objective algorithm. @see stop_condition @see StopCondition \n
        * Must be at least 1.
        *
        * @param count The number of generations to look back when checking the stall conditions.
        */
        void stall_gen_count(size_t count)
        {
            if (count == 0) throw std::invalid_argument("The stall generation count must be at least 1.");

            stall_gen_count_ = count;
        }
        [[nodiscard]] size_t stall_gen_count() const { return stall_gen_count_; }

        /**
        * Sets the value of the stall threshold to @p threshold for the stall stop conditions. \n
        * Only relevant for the single-objective algorithm. @see stop_condition @see StopCondition \n
        * May be negative if the deterioration of the stall metric can be allowed.
        *
        * @param threshold The stall threshold to use.
        */
        void stall_threshold(double threshold)
        {
            if (!std::isfinite(threshold)) throw std::invalid_argument("The stall threshold must be finite.");

            stall_threshold_ = threshold;
        }
        [[nodiscard]] double stall_threshold() const { return stall_threshold_; }

        /**
        * Sets the initial population to be used in the algorithm to @p pop instead of randomly generating it. \n
        * If @p pop is empty, the initial population will be randomly generated. \n
        * If the preset population's size is not equal to the population size set, either additional randomly generated
        * Candidates will be added to fill out the initial population, or some Candidates will be discarded from the preset.
        *
        * @param pop The initial population to use in the algorithm.
        */
        void presetInitialPopulation(const Population& pop)
        {
            if (!std::all_of(pop.begin(), pop.end(), [this](const Candidate& c) { return c.chromosome.size() == chrom_len_; }))
            {
                throw std::invalid_argument("The length of each chromosome in the preset pop must be equal to chrom_len.");
            }

            initial_population_preset_ = pop;
        }

        /**
        * Sets the fitness function used by the algorithm to @p f. \n
        * The fitness function should return a vector whose size is equal to the number of objectives, and
        * each element of the vector should be finite.
        *
        * @param f The fitness function to find the maximum of.
        */
        void setFitnessFunction(fitnessFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The fitness function can't be a nullptr.");

            fitnessFunction = f;
        }

        /* Some getters for the NSGA-III algorithm. */
        [[nodiscard]] std::vector<std::vector<double>> ref_points() const { return ref_points_; };
        [[nodiscard]] std::vector<double> ideal_point() const { return ideal_point_; };
        [[nodiscard]] std::vector<double> nadir_point() const { return nadir_point_; };

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

        void init()
        {
            /* Checks. */
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

        virtual Candidate generateCandidate() const = 0;

        Population generateInitialPopulation() const
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

        void evaluate(Population& pop)
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

        void updateOptimalSolutions(CandidateVec& optimal_sols, const Population& pop) const
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

        void prepSelections(Population& pop) const
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

        Candidate select(const Population& pop) const
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

        virtual CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const = 0;

        virtual void mutate(Candidate& child) const = 0;

        void repair(Population& pop) const
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

        Population updatePopulation(Population& old_pop, CandidateVec& children)
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

        bool stopCondition() const
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

        void updateStats(const Population& pop)
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

        /* SOGA functions. */

        /* Functions for calculating the selection probabilities of individuals in the single-objective algorithm. */

        static void sogaCalcRouletteWeights(Population& pop)
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

        static void sogaCalcRankWeights(Population& pop, double weight_min = 0.1, double weight_max = 1.1)
        {
            assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
            assert(0.0 <= weight_min && weight_min < weight_max && weight_max <= std::numeric_limits<double>::max());

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

        static void sogaCalcSigmaWeights(Population& pop, double scale = 3.0)
        {
            assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
            assert(scale > 1.0);

            double fitness_mean = fitnessMean(pop);
            double fitness_sd = fitnessSD(pop);

            double pdf_mean = 0.0;
            for (auto& sol : pop)
            {
                sol.selection_pdf = 1.0 + (sol.fitness[0] - fitness_mean) / (scale * std::max(fitness_sd, 1E-6));
                sol.selection_pdf = std::max(sol.selection_pdf, 0.0);	/* If (fitness < f_mean - scale * SD) the weight could be negative. */

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

        static void sogaCalcBoltzmannWeights(Population& pop, size_t t, size_t t_max, double temp_min, double temp_max)
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

        void sogaCalcWeights(Population& pop) const
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

        /* Functions used for the selections in the single-objective algorithm. */

        static Candidate sogaWeightProportionalSelect(const Population& pop)
        {
            assert(!pop.empty());
            assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

            double threshold = rng::generateRandomDouble();
            auto it = std::lower_bound(pop.begin(), pop.end(), threshold,
            [](const Candidate& sol, double threshold)
            {
                return sol.selection_cdf < threshold;
            });

            return (it != pop.end()) ? *it : pop.back();
        }

        static Candidate sogaTournamentSelect(const Population& pop, size_t tourney_size)
        {	
            assert(!pop.empty());
            assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));
            assert(tourney_size > 1);

            /* Randomly pick tourney_size candidates. Indices may repeat. */
            std::vector<size_t> indices;
            indices.reserve(tourney_size);
            for (size_t i = 0; i < tourney_size; i++)
            {
                indices.push_back(rng::generateRandomIdx(pop.size()));
            }

            /* Find the best of the picked candidates. */
            size_t idx = *std::max_element(indices.begin(), indices.end(),
            [&pop](size_t lidx, size_t ridx)
            {
                return pop[lidx].fitness[0] < pop[ridx].fitness[0];
            });

            return pop[idx];
        }

        Candidate sogaSelect(const Population& pop) const
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

        Population updateSogaPopulation(Population& old_pop, CandidateVec& children) const
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

        /* NSGA-II functions. */

        /* Find all Pareto fronts in the population and also assign the nondomination ranks of the candidates. */
        static std::vector<std::vector<size_t>> nonDominatedSort(Population& pop)
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

        /* Calculate the crowding distances of the candidates in each pareto front of the population. */
        static void calcCrowdingDistances(Population& pop, std::vector<std::vector<size_t>>& pfronts)
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

        /* Returns true if lhs is better than rhs. */
        static bool crowdedCompare(const Candidate& lhs, const Candidate& rhs)
        {
            if (rhs.rank > lhs.rank) return true;
            else if (lhs.rank == rhs.rank) return lhs.distance > rhs.distance;
            else return false;
        }

        /* Tournament selection using the crowding distances for tiebreaks. */
        static Candidate nsga2Select(const Population& pop)
        {
            assert(!pop.empty());

            size_t idx1 = rng::generateRandomIdx(pop.size());
            size_t idx2 = rng::generateRandomIdx(pop.size());

            return crowdedCompare(pop[idx1], pop[idx2]) ? pop[idx1] : pop[idx2];
        }

        Population updateNsga2Population(Population& old_pop, CandidateVec& children) const
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

        /* NSGA-III functions. */

        void updateIdealPoint(const Population& pop)
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

        void updateNadirPoint(const Population& pop)
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

        /* Find the closest reference point to each candidate after normalization, and their distances. */
        void associatePopToRefs(Population& pop, const std::vector<std::vector<double>>& ref_points)
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

        /* Return the niche counts of the ref points and assign niche counts to the candidates. */
        static std::vector<size_t> calcNicheCounts(Population& pop, const std::vector<std::vector<double>>& ref_points)
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

        /* Returns true if lhs is better than rhs. */
        static bool nichedCompare(const Candidate& lhs, const Candidate& rhs)
        {
            if (rhs.rank > lhs.rank) return true;
            else if (lhs.rank == rhs.rank) return lhs.niche_count < rhs.niche_count;
            else if (lhs.rank == rhs.rank && lhs.niche_count == rhs.niche_count) return lhs.distance < rhs.distance;
            else return false;
        }

        /* Tournament selection using the niche counts for tiebreaks. */
        static Candidate nsga3Select(const Population& pop)
        {
            assert(!pop.empty());

            size_t idx1 = rng::generateRandomIdx(pop.size());
            size_t idx2 = rng::generateRandomIdx(pop.size());

            return nichedCompare(pop[idx1], pop[idx2]) ? pop[idx1] : pop[idx2];
        }

        Population updateNsga3Population(Population& old_pop, CandidateVec& children)
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
                size_t ref = refs[rng::generateRandomIdx(refs.size())];

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

        /* Functions used to find the optimal solutions in the population. */

        static CandidateVec findParetoFront1D(const Population& pop)
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

        static CandidateVec findParetoFrontKung(const Population& pop)
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

        /* Utility functions. */

        /* Find the minimum fitness along each objective in the population. */
        static std::vector<double> fitnessMin(const Population& pop)
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

        /* Find the maximum fitness along each objective in the population. */
        static std::vector<double> fitnessMax(const Population& pop)
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

        /* Find the mean of the fitness values of the population along the first objective. */
        static double fitnessMean(const Population& pop)
        {
            assert(!pop.empty());
            assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return !sol.fitness.empty(); }));

            return std::accumulate(pop.begin(), pop.end(), 0.0,
            [&pop](double sum, const Candidate& sol)
            {
                return sum + sol.fitness[0] / pop.size(); 
            });
        }

        /* Find the standard deviation of the fitness values of the population along the first objective. */
        static double fitnessSD(const Population& pop)
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

    };

    /**
    * Standard genetic algorithm with binary encoding. \n
    * (The binary genes are encoded as char types.)
    */
    class BinaryGA : public GA<char>
    {
    public:

        /**
        * Possible crossover operators that can be used in the BinaryGA. \n
        * Set the crossover method used in the algorithm with @ref crossover_method. \n
        * The function used for the crossovers with the custom method can be set with @ref setCrossoverFunction.
        */
        enum class CrossoverMethod
        {
            single_point,	/**< Single-point crossover operator. */
            two_point,		/**< Two-point crossover operator. */
            n_point,		/**< General n-point crossover operator. @see num_crossover_points. */
            uniform,		/**< Uniform crossover operator. */
            custom			/**< Custom crossover operator defined by the user. @see setCrossoverFunction */
        };

        /**
        * Possible mutation operators that can be used in the BinaryGA. \n
        * Set the mutation method used in the algorithm with @ref mutation_method. \n
        * The function used for the mutations with the custom method can be set with @ref setMutationFunction.
        */
        enum class MutationMethod
        {
            standard,		/**< Standard mutation operator used in binary coded genetic algorithms. */
            custom			/**< Custom mutation operator defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the binary GA.
        * 
        * @param chrom_len The length of the binary chromosomes.
        * @param fitness_function The fitness function to find the maximum of in the algorithm.
        */
        BinaryGA(size_t chrom_len, fitnessFunction_t fitness_function) : GA(chrom_len, fitness_function) {}

        /**
        * Sets the crossover function used in the algorithm to @f.
        * @see CrossoverMethod
        *
        * @param method The crossover function to use.
        */
        void crossover_method(crossoverFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            crossover_method_ = CrossoverMethod::custom;
            customCrossover = f;
        }

        /**
        * Sets the crossover method used in the algorithm to @p method.
        * @see CrossoverMethod
        *
        * @param method The crossover method to use.
        */
        void crossover_method(CrossoverMethod method)
        {
            if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

            crossover_method_ = method;
        }
        [[nodiscard]] CrossoverMethod crossover_method() const { return crossover_method_; }

        /**
        * Sets the mutation function used in the algorithm to @f.
        * @see MutationMethod
        *
        * @param method The mutation function to use.
        */
        void mutation_method(mutationFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            mutation_method_ = MutationMethod::custom;
            customMutate = f;
        }

        /**
        * Sets the mutation method used in the algorithm to @p method.
        * @see MutationMethod
        *
        * @param method The mutation method to use.
        */
        void mutation_method(MutationMethod method)
        {
            if (static_cast<size_t>(method) > 1) throw std::invalid_argument("Invalid mutation method selected.");

            mutation_method_ = method;
        }
        [[nodiscard]] MutationMethod mutation_method() const { return mutation_method_; }

        /**
        * Sets the number of crossover points used in the crossovers to @p n if the n_point crossover method selected. \n
        * The number of crossover points must be at least 1.
        * @see crossover_method @see CrossoverMethod
        *
        * @param n The number of crossover points.
        */
        void num_crossover_points(size_t n)
        {
            if (n == 0) throw std::invalid_argument("The number of crossover points must be at least 1.");

            num_crossover_points_ = n;
        }
        [[nodiscard]] size_t num_crossover_points() const { return num_crossover_points_; }

    private:
        
        /* Parameters specific to the binary GA. */

        CrossoverMethod crossover_method_ = CrossoverMethod::single_point;
        MutationMethod mutation_method_ = MutationMethod::standard;
        size_t num_crossover_points_ = 3;

        Candidate generateCandidate() const override
        {
            assert(chrom_len_ > 0);

            Candidate sol;
            sol.chromosome.reserve(chrom_len_);
            for (size_t i = 0; i < chrom_len_; i++)
            {
                sol.chromosome.push_back(char(rng::generateRandomBool()));
            }

            return sol;
        }

        /* Funtions used for the crossovers. */

        static CandidatePair nPointCrossover(const Candidate& parent1, const Candidate& parent2, double pc, size_t n)
        {	
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Generate n (or less, but at least 1) number of unique random loci. */
                unordered_set<size_t> loci;
                for (size_t i = 0; i < n; i++)
                {
                    loci.insert(rng::generateRandomInt(size_t{ 1 }, parent1.chromosome.size() - 1));
                }

                /* Count how many loci are after each gene. */
                vector<size_t> loci_after;
                loci_after.reserve(parent1.chromosome.size());

                size_t loci_left = loci.size();
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    if (loci_left > 0 && loci.contains(i)) loci_left--;
                    loci_after.push_back(loci_left);
                }

                /* Perform the crossover. */
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* Swap the bits if there are an odd number of loci after the gene. */
                    if (loci_after[i] % 2)
                    {
                        child1.chromosome[i] = parent2.chromosome[i];
                        child2.chromosome[i] = parent1.chromosome[i];
                    }
                }

                /* Check if the children will need evaluation. */
                if (child1 != parent1)
                {
                    child1.is_evaluated = false;
                    child2.is_evaluated = false;
                }
            }

            return make_pair(child1, child2);
        }

        static CandidatePair uniformCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* Swap each gene with 0.5 probability. */
                    if (rng::generateRandomBool())
                    {
                        child1.chromosome[i] = parent2.chromosome[i];
                        child2.chromosome[i] = parent1.chromosome[i];
                    }
                }
                /* Check if the children will need evaluation. */
                if (child1 != parent1)
                {
                    child1.is_evaluated = false;
                    child2.is_evaluated = false;
                }
            }

            return std::make_pair(child1, child2);
        }

        CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const override
        {
            /* Edge case. No point in performing the crossover if the parents are the same. */
            if (parent1 == parent2) return std::make_pair(parent1, parent2);

            switch (crossover_method_)				
            {
                case CrossoverMethod::single_point:
                    return nPointCrossover(parent1, parent2, crossover_rate_, 1);
                case CrossoverMethod::two_point:
                    return nPointCrossover(parent1, parent2, crossover_rate_, 2);
                case CrossoverMethod::n_point:
                    return nPointCrossover(parent1, parent2, crossover_rate_, num_crossover_points_);
                case CrossoverMethod::uniform:
                    return uniformCrossover(parent1, parent2, crossover_rate_);
                case CrossoverMethod::custom:
                    return customCrossover(parent1, parent2, crossover_rate_);
                default:
                    assert(false);	/* Invalid crossover method. Shouldn't get here. */
                    std::abort();
            }
        }

        /* Functions used for the mutations. */

        static void standardMutate(Candidate& child, double pm)
        {
            assert(0.0 <= pm && pm <= 1.0);

            /* Calc number of mutated genes. */
            double mean = child.chromosome.size() * pm;
            double SD = child.chromosome.size() * pm * (1.0 - pm);

            size_t mutation_count = size_t(round(rng::generateRandomNorm(mean, SD)));
            mutation_count = std::clamp(mutation_count, size_t{ 0 }, child.chromosome.size());
                
            /* The child will (very likely) be changed, and will need to be evaluated. */
            if (mutation_count > 0) child.is_evaluated = false;

            /* Flip mutation_count number of random genes. One gene may be flipped multiple times, but uncommon for long chromosomes. */
            for (size_t i = 0; i < mutation_count; i++)
            {
                size_t idx = rng::generateRandomIdx(child.chromosome.size());
                child.chromosome[idx] = !child.chromosome[idx];
            }
        }

        void mutate(Candidate& child) const override
        {
            switch (mutation_method_)
            {
                case MutationMethod::standard:
                    standardMutate(child, mutation_rate_);
                    break;
                case MutationMethod::custom:
                    customMutate(child, mutation_rate_);
                    break;
                default:
                    assert(false);	/* Invalid mutation method. Shouldnt get here. */
                    std::abort();
            }
        }
    };

    /**
    * Standard genetic algorithm that uses real encoding. \n
    * Each gene of the chromosomes is a real value.
    */
    class RCGA : public GA<double>
    {
    public:

        /**
        * For the gene boundaries. \n
        * For example: { {gene1_min, gene1_max}, {gene2_min, gene2_max}, ... }
        */
        using limits_t = std::vector<std::pair<double, double>>;

        /**
        * Possible crossover operators that can be used in the RCGA. \n
        * Includes commonly used crossover operators in real-coded genetic algorithms, but a custom function
        * can also be used to perform the crossovers with the custom option. \n
        * Set the crossover method used in the algorithm with @ref crossover_method.
        */
        enum class CrossoverMethod
        {
            arithmetic,			/**< Arithmetic crossover operator.	Uses no parameters. */
            blx_a,				/**< BLX-alpha (blend) crossover operator. @see blx_crossover_param */
            simulated_binary,	/**< Simulated binary crossover (SBX) operator. @see sim_binary_crossover_param */
            wright,				/**< Wright heuristic crossover (HX) operator. Uses no parameters. */
            custom				/**< Custom crossover operator defined by the user. @see setCrossoverFunction */
        };

        /**
        * Possible mutation operators that can be used in the RCGA. \n
        * Includes commonly used mutation operators in real-coded genetic algorithms, but a custom function
        * can also be used to perform the mutations with the custom option. \n
        * Set the mutation method used in the algorithm with @ref mutation_method.
        */
        enum class MutationMethod
        {
            random,			/**<  Random (uniform) mutation operator. Uses no parameters. */
            polynomial,		/**<  Polynomial mutation operator. @see polynomial_mutation_param */
            nonuniform,		/**<  Non-uniform mutation operator. @see nonuniform_mutation_param */
            boundary,		/**<  Boundary mutation operator. Uses no parameters. */
            gauss,			/**<  Gauss mutation operator. @see gauss_mutation_param */
            custom			/**<  Custom mutation operator defined by the user. Uses the @ref customMutate to perform the mutations. */
        };

        /**
        * Basic constructor for the RCGA.
        * 
        * @param chrom_len The number of real genes in the chromosomes of the candidates.
        * @param fitness_function The fitness function to find the maximum of with the algorithm.
        * @param bounds The boundaries of the real coded genes (their min and max values).
        */
        RCGA(size_t chrom_len, fitnessFunction_t fitnessFunction, limits_t bounds) : GA(chrom_len, fitnessFunction), limits_(bounds)
        {
            if (bounds.size() != chrom_len)
            {
                throw std::invalid_argument("The size of the bounds must be the same as the number of genes.");
            }
            if (std::any_of(bounds.begin(), bounds.end(), [](std::pair<double, double> b) { return b.first > b.second; }))
            {
                throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
            }

        }

        /** Comparison operator for the real-coded candidates. */
        friend bool operator==(const Candidate& lhs, const Candidate& rhs)
        {
            return std::equal(lhs.chromosome.begin(), lhs.chromosome.end(), rhs.chromosome.begin(),
            [](double lhs, double rhs) 
            {
                return std::abs(lhs - rhs) <= std::numeric_limits<double>::epsilon() * std::max(std::abs(lhs), std::abs(rhs));
            });
        }

        /**
        * Sets the crossover function used in the algorithm to @f.
        * @see CrossoverMethod
        *
        * @param method The crossover function to use.
        */
        void crossover_method(crossoverFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            crossover_method_ = CrossoverMethod::custom;
            customCrossover = f;
        }

        /**
        * Sets the crossover method used in the algorithm to @p method.
        * @see CrossoverMethod
        *
        * @param method The crossover method to use.
        */
        void crossover_method(CrossoverMethod method)
        {
            if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

            crossover_method_ = method;
        }
        [[nodiscard]] CrossoverMethod crossover_method() const { return crossover_method_; }

        /**
        * Sets the mutation function used in the algorithm to @f.
        * @see MutationMethod
        *
        * @param method The mutation function to use.
        */
        void mutation_method(mutationFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            mutation_method_ = MutationMethod::custom;
            customMutate = f;
        }

        /**
        * Sets the mutation method used in the algorithm to @p method.
        * @see MutationMethod
        *
        * @param method The mutation method to use.
        */
        void mutation_method(MutationMethod method)
        {
            if (static_cast<size_t>(method) > 5) throw std::invalid_argument("Invalid mutation method selected.");

            mutation_method_ = method;
        }
        [[nodiscard]] MutationMethod mutation_method() const { return mutation_method_; }

        /**
        * Sets the boundaries of the real-coded genes. \n
        * Each element must contain the lower and upper bounds of the corresponding real gene. (The min and max values of the gene.) \n
        * The number of elements must be the same as the length of the chromosomes, and the lower bounds must not be higher than the upper bounds. \n
        * Eg. in the case of chromosomes of 2 values where both genes must be between -1 and 1: \n
        * limits = { {-1.0, 1.0}, {-1.0, 1.0} }
        *
        * @param limits The boundaries of the real genes.
        */
        void limits(limits_t limits)
        {
            if (limits.size() != chrom_len_)
            {
                throw std::invalid_argument("The number of limits must be equal to the chromosome length.");
            }
            if (std::any_of(limits.begin(), limits.end(), [](std::pair<double, double> b) {return b.first > b.second; }))
            {
                throw std::invalid_argument("The lower bound must be lower than the upper bound for each gene.");
            }

            limits_ = limits;
        }
        [[nodiscard]] limits_t limits() const { return limits_; }

        /**
        * Sets the alpha parameter of the BLX-alpha crossover method. \n
        * The parameter controls the length of the interval the children's genes are randomly chosen from, larger alpha -> larger interval. \n
        * Must be nonnegative, with the ideal being around 0.5.
        * @see crossover_method @see CrossoverMethod
        *
        * @param alpha The alpha parameter of the BLX-alpha crossover.
        */
        void blx_crossover_param(double alpha)
        {
            if (!(0.0 <= alpha && alpha <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Alpha must be a nonnegative, finite value.");
            }

            blx_crossover_param_ = alpha;
        }
        [[nodiscard]] double blx_crossover_param() const { return blx_crossover_param_; }

        /**
        * Sets the shape parameter (eta) of the simulated binary crossover method. \n
        * Must be nonnegative, typical values are 1-5. \n
        * Larger values will lead to the children's genes to be closer to the parents' genes.
        * @see crossover_method @see CrossoverMethod
        *
        * @param eta The parameter of the simulated binary crossover.
        */
        void sim_binary_crossover_param(double eta)
        {
            if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Eta must be a nonnegative, finite value.");
            }

            sim_binary_crossover_param_ = eta;
        }
        [[nodiscard]] double sim_binary_crossover_param() const { return sim_binary_crossover_param_; }

        /**
        * Sets the time parameter of the non-uniform mutation operator. \n
        * The mutated genes will be closer to the original values as the generations advance. Larger parameter values
        * will result in this process happening faster. \n
        * A value of 0 means the behaviour of the mutation does not change over time. \n
        * Must be nonnegative.
        * @see mutation_method @see MutationMethod
        *
        * @param b The parameter of the non-uniform mutation.
        */
        void nonuniform_mutation_param(double b)
        {
            if (!(0.0 <= b && b <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("The parameter b must be a nonnegative, finite value.");
            }

            nonuniform_mutation_param_ = b;
        }
        [[nodiscard]] double nonuniform_mutation_param() const { return nonuniform_mutation_param_; }

        /**
        * Sets the parameter of the polynomial mutation operator. \n
        * Must be nonnegative, typically in the range [20, 100]. \n
        * The mutated genes will be closer to their original values with higher parameter values.
        * @see mutation_method @see MutationMethod
        *
        * @param eta The parameter of the polynomial mutation.
        */
        void polynomial_mutation_param(double eta)
        {
            if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Eta must be a nonnegative, finite value.");
            }

            polynomial_mutation_param_ = eta;
        }
        [[nodiscard]] double polynomial_mutation_param() const { return polynomial_mutation_param_; }

        /**
        * Sets the parameter of the gauss mutation operator. \n
        * Must be larger than 0. \n
        * Controls how many standard deviations the gene's allowed range should take up. Larger values will result in mutated genes
        * closer to the original genes.
        * @see mutation_method @see MutationMethod
        *
        * @param sigmas The parameter of the gauss mutation.
        */
        void gauss_mutation_param(double sigmas)
        {
            if (!(0.0 < sigmas && sigmas <= std::numeric_limits<double>::max()))
            {
                throw std::invalid_argument("Eta must be a nonnegative, finite value.");
            }

            gauss_mutation_param_ = sigmas;
        }
        [[nodiscard]] double gauss_mutation_param() const { return gauss_mutation_param_; }

    private:

        /* Parameters specific to the real-coded GA. */

        limits_t limits_;

        CrossoverMethod crossover_method_ = CrossoverMethod::blx_a;
        double blx_crossover_param_ = 0.5;
        double sim_binary_crossover_param_ = 4.0;

        MutationMethod mutation_method_ = MutationMethod::random;
        double nonuniform_mutation_param_ = 2.0;
        double polynomial_mutation_param_ = 40.0;
        double gauss_mutation_param_ = 6.0;

        Candidate generateCandidate() const override
        {
            assert(chrom_len_ > 0);
            assert(chrom_len_ == limits_.size());

            Candidate sol;
            sol.chromosome.reserve(chrom_len_);
            for (size_t i = 0; i < chrom_len_; i++)
            {
                sol.chromosome.push_back(rng::generateRandomDouble(limits_[i].first, limits_[i].second));
            }

            return sol;
        }

        /* Function used for the crossovers. */

        static CandidatePair arithmeticCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                double alpha = rng::generateRandomDouble();
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    child1.chromosome[i] = alpha * parent1.chromosome[i] + (1.0 - alpha) * parent2.chromosome[i];
                    child2.chromosome[i] = (1.0 - alpha) * parent1.chromosome[i] + alpha * parent2.chromosome[i];
                }
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            return std::make_pair(child1, child2);
        }

        static CandidatePair blxAlphaCrossover(const Candidate& parent1, const Candidate& parent2, double pc, double alpha, const limits_t& bounds)
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(parent1.chromosome.size() == bounds.size());
            assert(0.0 <= pc && pc <= 1.0);
            assert(alpha >= 0.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* Calc interval to generate the childrens genes on. */
                    auto [range_min, range_max] = std::minmax(parent1.chromosome[i], parent2.chromosome[i]);
                    double range_ext = alpha * (range_max - range_min);
                    /* Generate genes from an uniform distribution on the interval. */
                    child1.chromosome[i] = rng::generateRandomDouble(range_min - range_ext, range_max + range_ext);
                    child2.chromosome[i] = rng::generateRandomDouble(range_min - range_ext, range_max + range_ext);
                    /* The generated genes might be outside the allowed interval. */
                    child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].first, bounds[i].second);
                    child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].first, bounds[i].second);
                }
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            return std::make_pair(child1, child2);
        }

        static CandidatePair simulatedBinaryCrossover(const Candidate& parent1, const Candidate& parent2, double pc, double b, const limits_t& bounds)
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(parent1.chromosome.size() == bounds.size());
            assert(0.0 <= pc && pc <= 1.0);
            assert(b > 0.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Generate beta from the distribution. */
                double u = rng::generateRandomDouble(0.0, 1.0);
                double beta = (u <= 0.5) ? std::pow(2.0 * u, 1.0/(b + 1.0)) : std::pow(1.0/(2.0 * (1.0 - u)), 1.0/(b + 1.0));

                /* Perform crossover. */
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    child1.chromosome[i] = 0.5 * ((1 - beta) * parent1.chromosome[i] + (1 + beta) * parent2.chromosome[i]);
                    child2.chromosome[i] = 0.5 * ((1 + beta) * parent1.chromosome[i] + (1 - beta) * parent2.chromosome[i]);
                    /* The childrens genes may be outside the allowed interval. */
                    child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].first, bounds[i].second);
                    child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].first, bounds[i].second);
                }
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            return std::make_pair(child1, child2);
        }

        static CandidatePair wrightCrossover(const Candidate& parent1, const Candidate& parent2, double pc, const limits_t& bounds)
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(parent1.chromosome.size() == bounds.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* p1 is always the better parent. */
                const Candidate* p1 = detail::paretoCompare(parent1.fitness, parent2.fitness) ? &parent2 : &parent1;
                const Candidate* p2 = detail::paretoCompare(parent1.fitness, parent2.fitness) ? &parent1 : &parent2;
                /* Get random weights. */
                double w1 = rng::generateRandomDouble();
                double w2 = rng::generateRandomDouble();
                /* Perform crossover. */
                for (size_t i = 0; i < p1->chromosome.size(); i++)
                {
                    child1.chromosome[i] = w1 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
                    child2.chromosome[i] = w2 * (p1->chromosome[i] - p2->chromosome[i]) + p1->chromosome[i];
                    /* The generated childrens genes may be outside the allowed intervals. */
                    child1.chromosome[i] = std::clamp(child1.chromosome[i], bounds[i].first, bounds[i].second);
                    child2.chromosome[i] = std::clamp(child2.chromosome[i], bounds[i].first, bounds[i].second);
                }
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            return std::make_pair(child1, child2);
        }

        CandidatePair crossover(const Candidate& p1, const Candidate& p2) const override
        {
            switch (crossover_method_)
            {
                case CrossoverMethod::arithmetic:
                    return arithmeticCrossover(p1, p2, crossover_rate_);
                case CrossoverMethod::blx_a:
                    return blxAlphaCrossover(p1, p2, crossover_rate_, blx_crossover_param_, limits_);
                case CrossoverMethod::simulated_binary:
                    return simulatedBinaryCrossover(p1, p2, crossover_rate_, sim_binary_crossover_param_, limits_);
                case CrossoverMethod::wright:
                    return wrightCrossover(p1, p2, crossover_rate_, limits_);
                case CrossoverMethod::custom:
                    return customCrossover(p1, p2, crossover_rate_);
                default:
                    assert(false);	/* Invalid crossover method. Shouldn't get here. */
                    std::abort();
            }
        }

        /* Functions used for the mutations. */

        static void randomMutate(Candidate& child, double pm, const limits_t& bounds)
        {
            assert(0.0 <= pm && pm <= 1.0);
            assert(child.chromosome.size() == bounds.size());

            for (size_t i = 0; i < child.chromosome.size(); i++)
            {	
                /* Mutate the gene with pm probability. */
                if (rng::generateRandomDouble() <= pm)
                {
                    child.chromosome[i] = rng::generateRandomDouble(bounds[i].first, bounds[i].second);
                    child.is_evaluated = false;
                }
            }
        }

        static void nonuniformMutate(Candidate& child, double pm, size_t time, size_t time_max, double b, const limits_t& bounds)
        {
            assert(0.0 <= pm && pm <= 1.0);
            assert(child.chromosome.size() == bounds.size());
            assert(b >= 0.0);

            for (size_t i = 0; i < child.chromosome.size(); i++)
            {
                /* Perform mutation on the gene with pm probability. */
                if (rng::generateRandomDouble() <= pm)
                {
                    double interval = bounds[i].second - bounds[i].first;
                    double r = rng::generateRandomDouble();
                    double sign = rng::generateRandomBool() ? 1.0 : -1.0;

                    child.chromosome[i] += sign * interval * (1.0 - std::pow(r, std::pow(1.0 - double(time) / time_max, b)));
                    child.is_evaluated = false;

                    /* The mutated gene might be outside the allowed range. */
                    child.chromosome[i] = std::clamp(child.chromosome[i], bounds[i].first, bounds[i].second);
                }
            }
        }

        static void polynomialMutate(Candidate& child, double pm, double eta, const limits_t& bounds)
        {
            assert(0.0 <= pm && pm <= 1.0);
            assert(child.chromosome.size() == bounds.size());
            assert(eta >= 0.0);

            for (size_t i = 0; i < child.chromosome.size(); i++)
            {
                /* Perform mutation on the gene with pm probability. */
                if (rng::generateRandomDouble() <= pm)
                {
                    double u = rng::generateRandomDouble();
                    if (u <= 0.5)
                    {
                        double delta = std::pow(2.0 * u, 1.0 / (1.0 + eta)) - 1.0;
                        child.chromosome[i] += delta * (child.chromosome[i] - bounds[i].first);
                    }
                    else
                    {
                        double delta = 1.0 - std::pow(2.0 - 2.0 * u, 1.0 / (1.0 + eta));
                        child.chromosome[i] += delta * (bounds[i].second - child.chromosome[i]);
                    }
                    child.is_evaluated = false;
                    /* The mutated gene will always be in the allowed range. */
                }
            }
        }

        static void boundaryMutate(Candidate& child, double pm, const limits_t& bounds)
        {
            assert(0.0 <= pm && pm <= 1.0);
            assert(child.chromosome.size() == bounds.size());

            for (size_t i = 0; i < child.chromosome.size(); i++)
            {
                /* Perform mutation on the gene with pm probability. */
                if (rng::generateRandomDouble() <= pm)
                {
                    child.chromosome[i] = rng::generateRandomBool() ? bounds[i].first : bounds[i].second;
                    child.is_evaluated = false;
                }
            }
        }

        static void gaussMutate(Candidate& child, double pm, double scale, const limits_t& bounds)
        {
            assert(0.0 <= pm && pm <= 1.0);
            assert(child.chromosome.size() == bounds.size());
            assert(scale > 0.0);

            for (size_t i = 0; i < child.chromosome.size(); i++)
            {
                /* Perform mutation on the gene with pm probability. */
                if (rng::generateRandomDouble() <= pm)
                {
                    double SD = (bounds[i].second - bounds[i].first) / scale;
                    child.chromosome[i] += rng::generateRandomNorm(0.0, SD);
                    child.is_evaluated = false;
                    /* The mutated gene might be outside the allowed range. */
                    child.chromosome[i] = std::clamp(child.chromosome[i], bounds[i].first, bounds[i].second);
                }
            }
        }

        void mutate(Candidate& child) const override
        {
            switch (mutation_method_)
            {
                case MutationMethod::random:
                    randomMutate(child, mutation_rate_, limits_);
                    break;
                case MutationMethod::nonuniform:
                    nonuniformMutate(child, mutation_rate_, generation_cntr_, max_gen_, nonuniform_mutation_param_, limits_);
                    break;
                case MutationMethod::polynomial:
                    polynomialMutate(child, mutation_rate_, polynomial_mutation_param_, limits_);
                    break;
                case MutationMethod::boundary:
                    boundaryMutate(child, mutation_rate_, limits_);
                    break;
                case MutationMethod::gauss:
                    gaussMutate(child, mutation_rate_, gauss_mutation_param_, limits_);
                    break;
                case MutationMethod::custom:
                    customMutate(child, mutation_rate_);
                    break;
                default:
                    assert(false);	/* Invalid mutation method. Shouldnt get here. */
                    std::abort();
            }
        }
    };

    /**
    * Genetic algorithm that uses permutational encoding. \n
    * The genes of the chromosomes are all unique unsigned integers on [0, chrom_len-1].
    */
    class PermutationGA : public GA<size_t>
    {
    public:

        /**
        * Possible crossover methods that can be used in the PermutationGA. \n
        * Includes some commonly used crossover operators in permutations GAs, but a custom function can also be used
        * to perform the crossovers with the custom option. \n
        * Set the crossover method used in the algorithm with @ref crossover_method.
        */
        enum class CrossoverMethod
        {
            order,	/**< Order crossover operator (OX1). Uses no parameters. Fastest method. */
            cycle,	/**< Cycle crossover operator (CX). Uses no parameters. */
            edge,	/**< Edge assembly crossover operator (EAX). Uses no parameters. Slowest method. */
            pmx,	/**< Partially mapped crossover operator (PMX). Uses no parameters. */
            custom	/**< Custom crossover function defined by the user. @see setCrossoverFunction */
        };

        /**
        * Possible mutation methods that can be used in the PermutationGA. \n
        * Includes commonly used mutation operators in permutations GAs, but a custom mutation function can
        * also be used to perform the mutations with the custom option. \n
        * Set the mutation method used in the algorithm with @ref mutation_method.
        */
        enum class MutationMethod
        {
            swap,		/**< Single-swap mutation operator. Uses no parameters. */
            scramble,	/**< Scramble mutation operator. Uses no parameters. */
            inversion,	/**< Inversion mutation operator. Uses no parameters. */
            custom		/**< Custom mutation function defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the PermutationGA.
        * 
        * @param chrom_len The number of genes in the chromosomes.
        * @param fitness_function The fitness function used in the algorithm to find the maximum of.
        */
        PermutationGA(size_t chrom_len, fitnessFunction_t fitnessFunction) : GA(chrom_len, fitnessFunction) {}

        /**
        * Sets the crossover function used in the algorithm to @f.
        * @see CrossoverMethod
        *
        * @param method The crossover function to use.
        */
        void crossover_method(crossoverFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            crossover_method_ = CrossoverMethod::custom;
            customCrossover = f;
        }

        /**
        * Sets the crossover method used in the algorithm to @p method.
        * @see CrossoverMethod
        *
        * @param method The crossover method to use.
        */
        void crossover_method(CrossoverMethod method)
        {
            if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

            crossover_method_ = method;
        }
        [[nodiscard]] CrossoverMethod crossover_method() const { return crossover_method_; }

        /**
        * Sets the mutation function used in the algorithm to @f.
        * @see MutationMethod
        *
        * @param method The mutation function to use.
        */
        void mutation_method(mutationFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            mutation_method_ = MutationMethod::custom;
            customMutate = f;
        }

        /**
        * Sets the mutation method used in the algorithm to @p method.
        * @see MutationMethod
        *
        * @param method The mutation method to use.
        */
        void mutation_method(MutationMethod method)
        {
            if (static_cast<size_t>(method) > 3) throw std::invalid_argument("Invalid mutation method selected.");

            mutation_method_ = method;
        }
        [[nodiscard]] MutationMethod mutation_method() const { return mutation_method_; }

    private:

        /* PermutationGA specific parameters. */

        CrossoverMethod crossover_method_ = CrossoverMethod::order;
        MutationMethod mutation_method_ = MutationMethod::inversion;

        Candidate generateCandidate() const override
        {
            assert(chrom_len_ > 0);

            Candidate sol;
            static thread_local rng::PRNG engine{ std::random_device{}() };

            std::vector<size_t> chrom(chrom_len_);
            std::iota(chrom.begin(), chrom.end(), 0U);	
            std::shuffle(chrom.begin(), chrom.end(), engine);

            sol.chromosome = chrom;
            return sol;
        }

        /* Functions used for the crossovers. */

        static CandidatePair orderCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1, child2;
            child1.chromosome.reserve(parent1.chromosome.size());
            child2.chromosome.reserve(parent1.chromosome.size());

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Pick a random range of genes. */
                size_t r1 = rng::generateRandomIdx(parent1.chromosome.size());
                size_t r2 = rng::generateRandomIdx(parent1.chromosome.size());
                const auto [idx1, idx2] = minmax(r1, r2);

                /* Edge case. The entire chromosomes are swapped. */
                if (idx1 == 0 && idx2 == parent1.chromosome.size() - 1) return make_pair(parent2, parent1);

                /* The range that will go from parent1 -> child1. (Not using the constructor is intentional.) */
                unordered_set<size_t> range1;
                for (auto gene = parent1.chromosome.begin() + idx1; gene != parent1.chromosome.begin() + idx2 + 1; gene++) range1.insert(*gene);
                /* The range that will go from parent2 -> child2. */
                unordered_set<size_t> range2;
                for (auto gene = parent2.chromosome.begin() + idx1; gene != parent2.chromosome.begin() + idx2 + 1; gene++) range2.insert(*gene);

                /* Gather genes not in the range from the other parent. */
                vector<size_t> seg1;	/* Segment gathered from parent2 -> child1. */
                seg1.reserve(parent2.chromosome.size() - idx2 + idx1 - 1);
                vector<size_t> seg2;	/* Segment gathered from parent1 -> child2. */
                seg2.reserve(parent1.chromosome.size() - idx2 + idx1 - 1);
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* If this gene of parent2 is not in the range from parent1, add it to seg1. */
                    if (!range1.contains(parent2.chromosome[i])) seg1.push_back(parent2.chromosome[i]);
                    /* If this gene of parent1 is not in the range from parent2, add it to seg2. */
                    if (!range2.contains(parent1.chromosome[i])) seg2.push_back(parent1.chromosome[i]);
                }

                /* Construct the children. child1 = seg1 + p1_range, child2 = seg2 + p2_range */
                child1.chromosome.insert(child1.chromosome.end(), seg1.begin(), seg1.begin() + idx1);
                child1.chromosome.insert(child1.chromosome.end(), parent1.chromosome.begin() + idx1, parent1.chromosome.begin() + idx2 + 1);
                child1.chromosome.insert(child1.chromosome.end(), seg1.begin() + idx1, seg1.end());

                child2.chromosome.insert(child2.chromosome.end(), seg2.begin(), seg2.begin() + idx1);
                child2.chromosome.insert(child2.chromosome.end(), parent2.chromosome.begin() + idx1, parent2.chromosome.begin() + idx2 + 1);
                child2.chromosome.insert(child2.chromosome.end(), seg2.begin() + idx1, seg2.end());
            }
            else  /* No crossover. */
            {
                child1 = parent1;
                child2 = parent2;
            }

            return make_pair(child1, child2);
        }

        static CandidatePair cycleCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Identify all cycles. */
                vector<vector<size_t>> cycles;
                /* Copies of parent chromosomes so they can be changed without issues. */
                vector<size_t> chrom1 = parent1.chromosome;
                vector<size_t> chrom2 = parent2.chromosome;
                while (!chrom1.empty())
                {
                    /* Identify 1 cycle. */
                    vector<size_t> cycle;
                    /* Always start the cycle at chrom1[0]. */
                    size_t pos = 0;
                    size_t top_val = chrom1[pos];
                    cycle.push_back(top_val);

                    while (chrom2[pos] != chrom1[0])
                    {
                        /* Look for the bottom value at this pos (chrom2[pos]) in chrom1. This is the new pos and top value. */
                        pos = static_cast<size_t>(find(chrom1.begin(), chrom1.end(), chrom2[pos]) - chrom1.begin());
                        top_val = chrom1[pos];
                        /* Add the new top value to the cycle. */
                        cycle.push_back(top_val);
                        /* Keep going until the bottom value at this pos isn't the cycle start value. (cycle complete) */
                    }

                    /* Delete the values in this cycle from chrom1 and chrom2. (Without changing the order of the remaining genes.) */
                    for (size_t i = 0; i < cycle.size(); i++)
                    {
                        chrom1.erase(find(chrom1.begin(), chrom1.end(), cycle[i]));
                        chrom2.erase(find(chrom2.begin(), chrom2.end(), cycle[i]));
                    }
                    /* Add this cycle to the cycles. */
                    cycles.push_back(cycle);
                }

                /* Construct the children from the cycles. */
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* Find which cycle has the gene. */
                    size_t cycle_num = 0;
                    for (size_t j = 0; j < cycles.size(); j++)
                    {
                        if (find(cycles[j].begin(), cycles[j].end(), parent1.chromosome[i]) != cycles[j].end())
                        {
                            cycle_num = j + 1;
                            break;
                        }
                    }
                    /* Even cycle genes are swapped parent1->child2 and parent2->child1. */
                    if (cycle_num % 2 == 0)
                    {
                        child1.chromosome[i] = parent2.chromosome[i];
                        child2.chromosome[i] = parent1.chromosome[i];
                    }
                    /* Odd cycles were already handled when initializing the children. */
                }
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            return make_pair(child1, child2);
        }

        static CandidatePair edgeCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            using namespace std;
            using NList = vector<unordered_set<size_t>>;

            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1, child2;
            child1.chromosome.reserve(parent1.chromosome.size());
            child2.chromosome.reserve(parent2.chromosome.size());

            /* Crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Construct neighbour list based on parents. The first and last genes are not neighbours. */
                size_t len = parent1.chromosome.size();
                NList nl1(len);
                /* Neighbours of first and last genes. */
                nl1[parent1.chromosome.front()] = { parent1.chromosome[1] };
                nl1[parent1.chromosome.back()] = { parent1.chromosome[len - 2] };
                nl1[parent2.chromosome.front()].insert(parent2.chromosome[1]);
                nl1[parent2.chromosome.back()].insert(parent2.chromosome[len - 2]);
                /* Neighbours of all other genes. */
                for (size_t i = 1; i < len - 1; i++)
                {
                    nl1[parent1.chromosome[i]].insert(parent1.chromosome[i + 1]);
                    nl1[parent1.chromosome[i]].insert(parent1.chromosome[i - 1]);

                    nl1[parent2.chromosome[i]].insert(parent2.chromosome[i + 1]);
                    nl1[parent2.chromosome[i]].insert(parent2.chromosome[i - 1]);
                }
                NList nl2 = nl1;	/* Copy for child2. */

                /* Generate child1. */
                size_t X = parent1.chromosome[0];
                vector<size_t> not_in_child = parent1.chromosome;
                while (child1.chromosome.size() != parent1.chromosome.size())
                {
                    /* Append X to the child, and remove X from all neighbour lists. */
                    child1.chromosome.push_back(X);
                    not_in_child.erase(remove(not_in_child.begin(), not_in_child.end(), X), not_in_child.end());
                    for (auto& elem : nl1) elem.erase(X);

                    /* Determine next X that will be added to the child. */
                    if (child1.chromosome.size() != parent1.chromosome.size())
                    {
                        /* If X's neighbour list is empty, X = random node not already in child. */
                        if (nl1[X].empty())
                        {
                            X = not_in_child[rng::generateRandomIdx(not_in_child.size())];
                        }
                        else /* X's neighbour list is not empty, X = neighbour of X with fewest neighbours (random if tie). */
                        {
                            /* Find X's neighbour with fewest neighbours. */
                            size_t nb = *min_element(nl1[X].begin(), nl1[X].end(),
                            [&nl1](const size_t& lhs, const size_t& rhs)
                            {
                                return (nl1[lhs].size() < nl1[rhs].size());
                            });
                            size_t min_neighbour_count = nl1[nb].size();

                            /* Determine possible nodes (neighbours of X with min_neighbour_count num of neighbours). */
                            vector<size_t> possible_nodes;
                            for (const auto& elem : nl1[X])
                            {
                                if (nl1[elem].size() == min_neighbour_count) possible_nodes.push_back(elem);
                            }

                            X = possible_nodes[rng::generateRandomIdx(possible_nodes.size())];
                        }
                    }
                }
                /* Do same to get child2. */
                X = parent2.chromosome[0];
                not_in_child = parent2.chromosome;
                while (child2.chromosome.size() != parent2.chromosome.size())
                {
                    child2.chromosome.push_back(X);
                    not_in_child.erase(remove(not_in_child.begin(), not_in_child.end(), X), not_in_child.end());
                    for (auto& elem : nl2) { elem.erase(X); }

                    if (child2.chromosome.size() != parent2.chromosome.size())
                    {
                        if (nl2[X].empty())
                        {
                            X = not_in_child[rng::generateRandomIdx(not_in_child.size())];
                        }
                        else
                        {
                            size_t nb = *min_element(nl2[X].begin(), nl2[X].end(),
                            [&nl2](const size_t& lhs, const size_t& rhs)
                            {
                                return (nl2[lhs].size() < nl2[rhs].size());
                            });
                            size_t min_neighbour_count = nl2[nb].size();

                            vector<size_t> possible_nodes;
                            for (const auto& elem : nl2[X])
                            {
                                if (nl2[elem].size() == min_neighbour_count) possible_nodes.push_back(elem);
                            }
                            X = possible_nodes[rng::generateRandomIdx(possible_nodes.size())];
                        }
                    }
                }
            }
            else /* No crossover. */
            {
                child1 = parent1;
                child2 = parent2;
            }

            return make_pair(child1, child2);
        }

        static CandidatePair pmxCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);
            
            Candidate child1(parent2), child2(parent1);	/* Init so the last step of the crossover can be skipped. */

            /* Crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
                size_t r1 = rng::generateRandomIdx(parent1.chromosome.size());
                size_t r2 = rng::generateRandomIdx(parent1.chromosome.size());
                const auto [idx1, idx2] = minmax(r1, r2);

                /* Edge case. The entire chromosomes are copied directly. */
                if (idx1 == 0 && idx2 == parent1.chromosome.size() - 1) return make_pair(parent1, parent2);

                /* Copy values in the range from the corresponding parent. */
                for (size_t i = idx1; i <= idx2; i++)
                {
                    child1.chromosome[i] = parent1.chromosome[i];
                    child2.chromosome[i] = parent2.chromosome[i];
                }
                /* Ranges that were copied from parents for fast checking if they contain an element. (Not using the constructor is intentional.) */
                unordered_set<size_t> p1_range;
                for (auto gene = parent1.chromosome.begin() + idx1; gene != parent1.chromosome.begin() + idx2 + 1; gene++) p1_range.insert(*gene);
                unordered_set<size_t> p2_range;
                for (auto gene = parent2.chromosome.begin() + idx1; gene != parent2.chromosome.begin() + idx2 + 1; gene++) p2_range.insert(*gene);

                /* Get rest of the child genes from the other parents. */
                for (size_t i = idx1; i <= idx2; i++)
                {
                    /* Look for genes in parent2 in the same range which haven't been copied to child1 from parent1 (p1_range). */
                    if (!p1_range.contains(parent2.chromosome[i]))
                    {
                        size_t pos = i;
                        while (idx1 <= pos && pos <= idx2)
                        {
                            /* Look at value in parent1 in this same pos. */
                            size_t val = parent1.chromosome[pos];
                            /* Find this value in parent2. */
                            pos = static_cast<size_t>(find(parent2.chromosome.begin(), parent2.chromosome.end(), val) - parent2.chromosome.begin());
                            /* Keep going until pos is outside the range. */
                        }
                        child1.chromosome[pos] = parent2.chromosome[i];
                    }
                    /* Same for child2. */
                    if (!p2_range.contains(parent1.chromosome[i]))
                    {
                        size_t pos = i;
                        while (idx1 <= pos && pos <= idx2)
                        {
                            size_t val = parent2.chromosome[pos];
                            pos = static_cast<size_t>(find(parent1.chromosome.begin(), parent1.chromosome.end(), val) - parent1.chromosome.begin());
                        }
                        child2.chromosome[pos] = parent1.chromosome[i];
                    }
                }
                /* Copy any not yet in child positions to the children from the other parents. (Already done at the initialization of the children.) */
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            return make_pair(child1, child2);
        }

        CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const override
        {
            using namespace std;

            /* Edge case. No point in performing the crossover if the parents are the same. */
            if (parent1 == parent2)
            {
                return make_pair(parent1, parent2);
            }

            Candidate child1, child2;
            switch (crossover_method_)
            {
                case CrossoverMethod::order:
                    tie(child1, child2) = orderCrossover(parent1, parent2, crossover_rate_);
                    break;
                case CrossoverMethod::cycle:
                    tie(child1, child2) = cycleCrossover(parent1, parent2, crossover_rate_);
                    break;
                case CrossoverMethod::pmx:
                    tie(child1, child2) = pmxCrossover(parent1, parent2, crossover_rate_);
                    break;
                case CrossoverMethod::edge:
                    tie(child1, child2) = edgeCrossover(parent1, parent2, crossover_rate_);
                    break;
                case CrossoverMethod::custom:
                    tie(child1, child2) = customCrossover(parent1, parent2, crossover_rate_);
                    break;
                default:
                    assert(false);	/* Invalid crossover method. Shouldn't get here. */
                    abort();
            }

            /* Check if the evaluation of the children can be skipped. */
            /* 
            * These checks decrease fitness evals by a lot for short chromosomes:
            *	TSP13 (200pop, 1000gen, 0.9c):	~168'000 -> ~28'000 fitness evals -(second and last checks added)-> ~20'000 evals
            * Smaller decrease for long chromosomes:
            *	chrom_len=10'000 (50pop, 10gen, 1.0c): 500 -> ~480-495 fitness evals -(second and last checks added)-> ~475-490 evals
            */
            if (child1 == parent1)
            {
                child1.fitness = parent1.fitness;
                child1.is_evaluated = true;
            }
            else if (child1 == parent2)
            {
                child1.fitness = parent2.fitness;
                child1.is_evaluated = true;
            }
            if (child2 == parent2)
            {
                child2.fitness = parent2.fitness;
                child2.is_evaluated = true;
            }
            else if (child2 == parent1)
            {
                child2.fitness = parent1.fitness;
                child2.is_evaluated = true;
            }

            return make_pair(child1, child2);
        }

        /* Functions used for the mutations. */

        static void swapMutate(Candidate& child, double pm)
        {
            assert(0.0 <= pm && pm <= 1.0);

            /* Perform mutation with pm probability. */
            if (rng::generateRandomDouble() <= pm)
            {
                /* r1 and r2 might be the same index, but its rare for long chromosomes. */
                size_t r1 = rng::generateRandomIdx(child.chromosome.size());
                size_t r2 = rng::generateRandomIdx(child.chromosome.size());

                std::swap(child.chromosome[r1], child.chromosome[r2]);

                /* If the indices are different, the child was changed and will need evaluation. */
                if (r1 != r2) child.is_evaluated = false;
            }
        }

        static void scrambleMutate(Candidate& child, double pm)
        {
            assert(0.0 <= pm && pm <= 1.0);

            /* Perform mutation with pm probability. */
            if (rng::generateRandomDouble() <= pm)
            {
                /* Pick a random range of genes. The bounds may be the same, but its rare for long chromosomes. */
                size_t r1 = rng::generateRandomIdx(child.chromosome.size());
                size_t r2 = rng::generateRandomIdx(child.chromosome.size());
                auto [idx1, idx2] = std::minmax(r1, r2);

                static thread_local rng::PRNG engine{ std::random_device{}() };
                std::shuffle(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1, engine);

                /* If the indices are different, the child was very likely changed and will need evaluation. */
                if (r1 != r2) child.is_evaluated = false;
            }
        }

        static void inversionMutate(Candidate& child, double pm)
        {
            assert(0.0 <= pm && pm <= 1.0);

            /* Perform mutation with pm probability. */
            if (rng::generateRandomDouble() <= pm)
            {
                /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
                size_t r1 = rng::generateRandomIdx(child.chromosome.size());
                size_t r2 = rng::generateRandomIdx(child.chromosome.size());
                auto [idx1, idx2] = std::minmax(r1, r2);

                std::reverse(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1);

                /* If the indices are different, the child was changed and will need evaluation. */
                if (r1 != r2) child.is_evaluated = false;
            }
        }

        void mutate(Candidate& child) const override
        {
            switch (mutation_method_)
            {
                case MutationMethod::swap:
                    swapMutate(child, mutation_rate_);
                    break;
                case MutationMethod::scramble:
                    scrambleMutate(child, mutation_rate_);
                    break;
                case MutationMethod::inversion:
                    inversionMutate(child, mutation_rate_);
                    break;
                case MutationMethod::custom:
                    customMutate(child, mutation_rate_);
                    break;
                default:
                    assert(false);	/* Invalid mutation method. Shouldn't get here. */
                    std::abort();
            }
        }
    };

    /**
    * Integer coded GA. \n
    * Same as @ref BinaryGA, but the genes of the chromosomes can be any integer on [0, base], not just 0 or 1. \n
    * It also uses a slightly different mutation function with swaps and inversions.
    */
    class IntegerGA : public GA<size_t>
    {
    public:

        /**
        * Possible crossover operators that can be used in the IntegerGA. \n
        * These crossover method are the same as the ones used in the binary coded algorithm. \n
        * Set the crossover method used in the algorithm with @ref crossover_method. \n
        * The function used for the crossovers with the custom method can be set with @ref setCrossoverFunction.
        */
        enum class CrossoverMethod
        {
            single_point,	/**< Single-point crossover operator. */
            two_point,		/**< Two-point crossover operator. */
            n_point,		/**< General n-point crossover operator. Set the number of crossover points with @ref num_crossover_points */
            uniform,		/**< Uniform crossover operator. */
            custom			/**< Custom crossover operator defined by the user. @see setCrossoverMethod */
        };

        /**
        * Possible mutation operators that can be used in the IntegerGA. \n
        * Same operators as in the binary coded algorithm, with some small changes. \n
        * Set the mutation method used in the algorithm with @ref mutation_method. \n
        * The function used for the mutations with the custom method can be set with @ref setMutationFunction.
        */
        enum class MutationMethod
        {
            standard,	/**< Standard mutation operator used in the @ref BinaryGA with swap and inversion added. @see swap_rate @see inversion_rate */
            custom		/**< Custom mutation operator defined by the user. @see setMutationFunction */
        };

        /**
        * Basic contructor for the IntegerGA.
        * 
        * @param chrom_len The number of genes in each chromosome.
        * @param fitness_function The fitness function used in the algorithm.
        * @param base The number of values a gene can take. Must be > 1. If 2, same as the @ref BinaryGA.
        */
        IntegerGA(size_t chrom_len, fitnessFunction_t fitnessFunction, size_t base) : GA(chrom_len, fitnessFunction), base_(base) 
        {
            if (base < 2) throw std::invalid_argument("The base must be at least 2.");
        }

        /**
        * Sets the crossover function used in the algorithm to @f.
        * @see CrossoverMethod
        *
        * @param method The crossover function to use.
        */
        void crossover_method(crossoverFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            crossover_method_ = CrossoverMethod::custom;
            customCrossover = f;
        }

        /**
        * Sets the crossover method used in the algorithm to @p method.
        * @see CrossoverMethod
        *
        * @param method The crossover method to use.
        */
        void crossover_method(CrossoverMethod method)
        {
            if (static_cast<size_t>(method) > 4) throw std::invalid_argument("Invalid crossover method selected.");

            crossover_method_ = method;
        }
        [[nodiscard]] CrossoverMethod crossover_method() const { return crossover_method_; }

        /**
        * Sets the mutation function used in the algorithm to @f.
        * @see MutationMethod
        *
        * @param method The mutation function to use.
        */
        void mutation_method(mutationFunction_t f)
        {
            if (f == nullptr) throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

            mutation_method_ = MutationMethod::custom;
            customMutate = f;
        }

        /**
        * Sets the mutation method used in the algorithm to @p method.
        * @see MutationMethod
        *
        * @param method The mutation method to use.
        */
        void mutation_method(MutationMethod method)
        {
            if (static_cast<size_t>(method) > 1) throw std::invalid_argument("Invalid mutation method selected.");

            mutation_method_ = method;
        }
        [[nodiscard]] MutationMethod mutation_method() const { return mutation_method_; }

        /**
        * Sets the number of crossover points used in the crossovers to @p n if the n_point crossover method selected. \n
        * The number of crossover points must be at least 1.
        * @see crossover_method @see CrossoverMethod
        *
        * @param n The number of crossover points.
        */
        void num_crossover_points(size_t n)
        {
            if (n == 0) throw std::invalid_argument("The number of crossover points must be at least 1.");

            num_crossover_points_ = n;
        }
        [[nodiscard]] size_t num_crossover_points() const { return num_crossover_points_; }

        /**
        * Sets the number of values a gene can take to @p base. \n
        * The value of the base must be at least 2, and the GA is essentially the same as the
        * BinaryGA if the base is set to 2.
        * 
        * @param base The number of values a gene can be.
        */
        void base(size_t base)
        {
            if (base < 2) throw std::invalid_argument("The base must be at least 2.");

            base_ = base;
        }
        [[nodiscard]] size_t base() const { return base_; }

        /**
        * Sets the probability of a single swap occuring during the mutation of a Candidate to @p ps. \n
        * The value of ps must be on the closed interval [0.0, 1.0].
        * 
        * @param ps The probability of swap during mutation.
        */
        void swap_rate(double ps)
        {
            if (!(0.0 <= ps && ps <= 1.0)) throw std::invalid_argument("The probability of a swap must be in [0, 1].");

            swap_rate_ = ps;
        }
        [[nodiscard]] double swap_rate() const { return swap_rate_; }

        /**
        * Sets the probability of inversion during the mutation of a Candidate to @p pi. \n
        * The value of pi must be on the closed interval [0.0, 1.0].
        * 
        * @param pi The probability of inversion during mutation.
        */
        void inversion_rate(double pi)
        {
            if (!(0.0 <= pi && pi <= 1.0)) throw std::invalid_argument("The probability of inversion must be in [0, 1].");

            inversion_rate_ = pi;
        }
        [[nodiscard]] double inversion_rate() const { return inversion_rate_; }

    private:

        /* Parameters specific to the integer coded GA. */

        CrossoverMethod crossover_method_ = CrossoverMethod::single_point;
        MutationMethod mutation_method_ = MutationMethod::standard;
        size_t num_crossover_points_ = 3;
        size_t base_ = 4;
        double swap_rate_ = 0.1;
        double inversion_rate_ = 0.1;

        Candidate generateCandidate() const override
        {
            assert(chrom_len_ > 0);
            assert(base_ > 1);

            Candidate sol;
            sol.chromosome.reserve(chrom_len_);
            for (size_t i = 0; i < chrom_len_; i++)
            {
                sol.chromosome.push_back(rng::generateRandomInt(size_t{ 0 }, base_ - 1));
            }

            return sol;
        }

        /* Functions used for the crossovers. Same as the crossovers used in the binary GA. */

        static CandidatePair nPointCrossover(const Candidate& parent1, const Candidate& parent2, double pc, size_t n)
        {
            using namespace std;
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                /* Generate n (or less, but at least 1) number of random loci. */
                unordered_set<size_t> loci;
                for (size_t i = 0; i < n; i++)
                {
                    loci.insert(rng::generateRandomInt(size_t{ 1 }, parent1.chromosome.size() - 1));
                }
                
                /* Count how many loci are after each gene. */
                vector<size_t> loci_after;
                loci_after.reserve(parent1.chromosome.size());

                size_t loci_left = loci.size();
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    if ((loci_left > 0) && loci.contains(i)) loci_left--;
                    loci_after.push_back(loci_left);
                }

                /* Perform the crossover. */
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* Swap the bits if there are an odd number of loci after the gene. */
                    if (loci_after[i] % 2)
                    {
                        child1.chromosome[i] = parent2.chromosome[i];
                        child2.chromosome[i] = parent1.chromosome[i];
                    }
                }
                /* Check if the children will need evaluation. */
                if (child1 != parent1)
                {
                    child1.is_evaluated = false;
                    child2.is_evaluated = false;
                }
            }

            return make_pair(child1, child2);
        }

        static CandidatePair uniformCrossover(const Candidate& parent1, const Candidate& parent2, double pc)
        {
            assert(parent1.chromosome.size() == parent2.chromosome.size());
            assert(0.0 <= pc && pc <= 1.0);

            Candidate child1(parent1), child2(parent2);

            /* Perform crossover with pc probability. */
            if (rng::generateRandomDouble() <= pc)
            {
                for (size_t i = 0; i < parent1.chromosome.size(); i++)
                {
                    /* Swap each bit with 0.5 probability. */
                    if (rng::generateRandomBool())
                    {
                        child1.chromosome[i] = parent2.chromosome[i];
                        child2.chromosome[i] = parent1.chromosome[i];
                    }
                }
                /* Check if the children will need evaluation. */
                if (child1 != parent1)
                {
                    child1.is_evaluated = false;
                    child2.is_evaluated = false;
                }
            }

            return std::make_pair(child1, child2);
        }

        CandidatePair crossover(const Candidate& parent1, const Candidate& parent2) const override
        {
            /* Edge case. No point in performing the mutations if the parents are the same. */
            if (parent1 == parent2) return std::make_pair(parent1, parent2);

            switch (crossover_method_)
            {
                case CrossoverMethod::single_point:
                    return nPointCrossover(parent1, parent2, crossover_rate_, 1);
                case CrossoverMethod::two_point:
                    return nPointCrossover(parent1, parent2, crossover_rate_, 2);
                case CrossoverMethod::n_point:
                    return nPointCrossover(parent1, parent2, crossover_rate_, num_crossover_points_);
                case CrossoverMethod::uniform:
                    return uniformCrossover(parent1, parent2, crossover_rate_);
                case CrossoverMethod::custom:
                    return customCrossover(parent1, parent2, crossover_rate_);
                default:
                    assert(false);	/* Invalid crossover method. Shouldn't get here. */
                    std::abort();
            }
        }

        /* Functions used for the mutations. */

        static void standardMutate(Candidate& child, double pm, double ps, double pi, size_t base_)
        {
            assert(0.0 <= pm && pm <= 1.0);
            assert(0.0 <= ps && ps <= 1.0);
            assert(0.0 <= pi && pi <= 1.0);
            assert(base_ > 1);

            /* Calc number of mutated genes. */
            double mean = child.chromosome.size() * pm;
            double SD = child.chromosome.size() * pm * (1.0 - pm);

            size_t mutation_count = size_t(std::round(rng::generateRandomNorm(mean, SD)));
            mutation_count = std::clamp(mutation_count, size_t{ 0 }, child.chromosome.size());

            /* The child will (very likely) be changed, and will need to be evaluated. */
            if (mutation_count > 0) child.is_evaluated = false;

            /* 
            * Change mutation_count number of random genes.
            * One gene may be changed multiple times, but it's uncommon for long chromosomes, and not important.
            */
            for (size_t i = 0; i < mutation_count; i++)
            {
                size_t idx = rng::generateRandomIdx(child.chromosome.size());
                child.chromosome[idx] = rng::generateRandomInt(size_t{ 0 }, base_ - 1);
            }

            /* Perform swap with ps probability. */
            if (rng::generateRandomDouble() <= ps)
            {
                /* r1 and r2 might be the same index, but its rare for long chromosomes. */
                size_t r1 = rng::generateRandomIdx(child.chromosome.size());
                size_t r2 = rng::generateRandomIdx(child.chromosome.size());
                std::swap(child.chromosome[r1], child.chromosome[r2]);

                if (child.chromosome[r1] != child.chromosome[r2]) child.is_evaluated = false;
            }

            /* Perform inversion with pi probability. */
            if (rng::generateRandomDouble() <= pi)
            {
                /* Pick a random range of genes. The bounds of the range may be the same, but its rare for long chromosomes. */
                size_t r1 = rng::generateRandomIdx(child.chromosome.size());
                size_t r2 = rng::generateRandomIdx(child.chromosome.size());
                auto [idx1, idx2] = std::minmax(r1, r2);

                std::reverse(child.chromosome.begin() + idx1, child.chromosome.begin() + idx2 + 1);

                if (r1 != r2) child.is_evaluated = false;
            }
        }

        void mutate(Candidate& child) const override
        {
            switch (mutation_method_)
            {
                case MutationMethod::standard:
                    standardMutate(child, mutation_rate_, swap_rate_, inversion_rate_, base_);
                    break;
                case MutationMethod::custom:
                    customMutate(child, mutation_rate_);
                    break;
                default:
                    assert(false);	/* Invalid mutation method. Shouldn't get here. */
                    std::abort();
            }			
        }
    };

} // namespace genetic_algorithm

#endif // GENETIC_ALGORITHM_H