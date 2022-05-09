/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "single_objective.hpp"
#include "../algorithms/ga_info.hpp"
#include "../utility/rng.hpp"
#include <limits>
#include <utility>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::selection::single_objective
{
    void Roulette::prepare(const GaInfo&, const FitnessMatrix& pop)
    {
        auto selection_weights = dtl::rouletteWeights(pop);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    size_t Roulette::select(const GaInfo&, const FitnessMatrix&)
    {
        return rng::sampleCdf(cdf_);
    }

    Tournament::Tournament(size_t size)
    {
        this->size(size);
    }

    void Tournament::size(size_t size)
    {
        if (size < 2)
        {
            throw std::invalid_argument("The tournament size must be at least 2.");
        }

        tourney_size_ = size;
    }

    void Tournament::prepare(const GaInfo&, const FitnessMatrix&)
    { /* Nothing to do for tournament selection. */
    }

    size_t Tournament::select(const GaInfo&, const FitnessMatrix& pop)
    {
        assert(pop.size() >= tourney_size_);
        assert(std::all_of(pop.begin(), pop.end(), [](const FitnessVector& sol) { return sol.size() == 1; }));

        auto candidates = rng::sampleUnique(pop.size(), tourney_size_);

        return *std::max_element(candidates.begin(), candidates.end(),
        [&pop](size_t lidx, size_t ridx)
        {
            return pop[lidx][0] < pop[ridx][0];
        });
    }

    Rank::Rank(double min_weight, double max_weight)
    {
        this->weights(min_weight, max_weight);
    }

    void Rank::min_weight(double min_weight)
    {
        this->weights(min_weight, max_weight_);
    }

    void Rank::max_weight(double max_weight)
    {
        this->weights(min_weight_, max_weight);
    }

    void Rank::weights(double min_weight, double max_weight)
    {
        if (!(0.0 <= min_weight && min_weight <= max_weight))
        {
            throw std::invalid_argument("The minimum weight must be in the closed interval [0.0, max_weight].");
        }
        if (!(min_weight <= max_weight && max_weight <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("The maximum weight must be in the closed interval [min_weight, DBL_MAX].");
        }

        min_weight_ = min_weight;
        max_weight_ = max_weight;
    }

    void Rank::prepare(const GaInfo&, const FitnessMatrix& pop)
    {
        auto selection_weights = dtl::rankWeights(pop, min_weight_, max_weight_);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    size_t Rank::select(const GaInfo&, const FitnessMatrix&)
    {
        return rng::sampleCdf(cdf_);
    }

    Sigma::Sigma(double scale)
    {
        this->scale(scale);
    }

    void Sigma::scale(double scale)
    {
        if (!(1.0 <= scale && scale <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Scale must be in the closed interval [1.0, DBL_MAX].");
        }

        scale_ = scale;
    }

    void Sigma::prepare(const GaInfo&, const FitnessMatrix& pop)
    {
        auto selection_weights = dtl::sigmaWeights(pop, scale_);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    size_t Sigma::select(const GaInfo&, const FitnessMatrix&)
    {
        return rng::sampleCdf(cdf_);
    }

    Boltzmann::Boltzmann(TemperatureFunction f)
        : temperature_(std::move(f))
    {
    }

    void Boltzmann::prepare(const GaInfo& ga, const FitnessMatrix& pop)
    {
        double T = temperature_(ga.generation_cntr(), ga.max_gen());
        auto selection_weights = dtl::boltzmannWeights(pop, T);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    size_t Boltzmann::select(const GaInfo&, const FitnessMatrix&)
    {
        return rng::sampleCdf(cdf_);
    }
} // namespace genetic_algorithm::selection::single_objective