#ifndef GA_BASE_GA_IMPL_HPP
#define GA_BASE_GA_IMPL_HPP

#include "base_ga.decl.hpp"
#include "crossover/crossover_base.hpp"
#include "mutation/mutation_base.hpp"
#include "stop_condition/stop_condition_base.hpp"
#include "rng.hpp"
#include "reference_points.h"
#include "mo_detail.h"
#include "math.hpp"
#include "population.hpp"

#include <execution>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <cassert>
#include <cmath>

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

    template<typename T>
    inline GA<T>::GA(size_t chrom_len, fitnessFunction_t fitness_function)
        : chrom_len_(chrom_len), fitnessFunction(fitness_function)
    {
        if (chrom_len == 0)
        {
            throw std::invalid_argument("The chromosome length must be at least 1."); // TODO assign only after check
        }
        if (fitnessFunction == nullptr)
        {
            throw std::invalid_argument("The fitness function is a nullptr.");
        }
    }

    template<typename T>
    inline typename const GA<T>::CandidateVec& GA<T>::solutions() const
    {
        return solutions_;
    }

    template<typename T>
    inline size_t GA<T>::num_fitness_evals() const
    {
        return static_cast<size_t>(num_fitness_evals_);
    }

    template<typename T>
    inline size_t GA<T>::generation_cntr() const
    {
        return generation_cntr_;
    }

    template<typename T>
    inline typename const GA<T>::Population& GA<T>::population() const
    {
        return population_;
    }

    template<typename geneType>
    inline typename const GA<geneType>::History& GA<geneType>::soga_history() const
    {
        return soga_history_;
    }

    template<typename GeneType>
    [[nodiscard]] size_t GA<GeneType>::num_objectives() const
    {
        return num_objectives_;
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

    template<typename GeneType>
    template<typename CrossoverType>
    //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>> && std::copy_constructible<CrossoverType>
    void GA<GeneType>::crossover_method(const CrossoverType& f)
    {
        crossover_ = std::make_unique<CrossoverType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::crossover_method(std::unique_ptr<crossover::Crossover<GeneType>>&& f)
    {
        if (f == nullptr) throw std::invalid_argument("The crossover method can't be a nullptr.");

        crossover_ = std::move(f);
    }

    template<typename GeneType>
    template<typename CrossoverType>
    //requires std::derived_from<CrossoverType, crossover::Crossover<GeneType>>
    CrossoverType& GA<GeneType>::crossover_method() const
    {
        return dynamic_cast<CrossoverType&>(*crossover_);
    }

    template<typename GeneType>
    template<typename MutationType>
    //requires std::derived_from<MutationType, mutation::Mutation<GeneType>> && std::copy_constructible<MutationType>
    void GA<GeneType>::mutation_method(const MutationType& f)
    {
        mutation_ = std::make_unique<MutationType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::mutation_method(std::unique_ptr<mutation::Mutation<GeneType>>&& f)
    {
        if (f == nullptr) throw std::invalid_argument("The mutation method can't be a nullptr.");

        mutation_ = std::move(f);
    }

    template<typename GeneType>
    template<typename MutationType>
    //requires std::derived_from<MutationType, mutation::Crossover<GeneType>>
    MutationType& GA<GeneType>::mutation_method() const
    {
        return dynamic_cast<MutationType&>(*mutation_);
    }

    template<typename GeneType>
    template<typename StopType>
    //requires std::derived_from<StopType, stopping::StopCondition<GA<GeneType>>> && std::copy_constructible<StopType>
    void GA<GeneType>::stop_condition(const StopType& f)
    {
        stop_condition_ = std::make_unique<StopType>(f);
    }

    template<typename GeneType>
    void GA<GeneType>::stop_condition(std::unique_ptr<stopping::StopCondition<GeneType>>&& f)
    {
        stop_condition_ = std::move(f);
    }

    template<typename GeneType>
    template<typename StopType>
    //requires std::derived_from<StopType, stopping::StopCondition<GeneType>>
    StopType& GA<GeneType>::stop_condition() const
    {
        return dynamic_cast<StopType&>(*stop_condition_);
    }

    template<typename GeneType>
    bool GA<GeneType>::stopCondition()
    {
        if (stop_condition_ != nullptr)
        {
            return (generation_cntr_ >= (max_gen_ - 1)) || (*stop_condition_)(*this);
        }
        else
        {
            return (generation_cntr_ >= (max_gen_ - 1));
        }
    }


    template<typename geneType>
    inline typename const GA<geneType>::CandidateVec& GA<geneType>::run()
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
                p = (*crossover_)(*this, p.first, p.second);
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
                (*mutation_)(*this, c);
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
        for (const auto& sol : optimal_sols) unique_sols.insert(std::move(sol));    /* Insertion is faster than ctor. */
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
                assert(false);    /* Invalid mode. Shouldn't get here. */
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
                assert(false);    /* Invalid mode. Shouldn't get here. */
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
                assert(false);    /* Invalid mode, shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline void GA<geneType>::updateStats(const Population& pop)
    {
        switch (mode_)
        {
            case Mode::single_objective:
                soga_history_.add(populationFitnessMean(pop)[0], populationFitnessSD(pop)[0], populationFitnessMin(pop)[0], populationFitnessMax(pop)[0]);
                break;
            case Mode::multi_objective_sorting:
                break;
            case Mode::multi_objective_decomp:
                break;
            default:
                assert(false);    /* Invalid mode, shouldn't get here. */
                std::abort();
        }
    }

    template<typename geneType>
    inline void GA<geneType>::sogaCalcRouletteWeights(Population& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        /* Roulette selection wouldn't work for negative fitness values. */
        bool has_negative_fitness = std::any_of(pop.begin(), pop.end(), [](const Candidate& sol) { return sol.fitness[0] < 0.0; });
        double offset = populationFitnessMin(pop)[0] * has_negative_fitness;

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

        double fitness_mean = populationFitnessMean(pop)[0];
        double fitness_sd = populationFitnessSD(pop)[0];

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

        double fmax = populationFitnessMax(pop)[0];
        double fmin = populationFitnessMin(pop)[0];

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
                assert(false);    /* Invalid selection method. Shouldn't get here. */
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
                assert(false);    /* Invalid selection method. Shouldn't get here. */
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
                if (detail::paretoCompareLess(pop[j].fitness, pop[i].fitness))
                {
                    dom_count[j]++;
                    dom_list[i].push_back(j);
                }
                else if (detail::paretoCompareLess(pop[i].fitness, pop[j].fitness))
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
            vector<size_t> added_indices(population_size_ - new_pop.size());    /* For updating the crowding distances in this front. */
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

        vector<vector<double>> fnorms(pop.size(), vector<double>(pop[0].fitness.size(), 0.0));    /* Don't change the actual fitness values. */

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

} // namespace genetic_algorithm

#endif // !GA_BASE_GA_IMPL_HPP