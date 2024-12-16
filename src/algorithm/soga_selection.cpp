/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "soga_selection.hpp"
#include "../core/ga_info.hpp"
#include "../core/candidate.hpp"
#include "../core/population.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/math.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <span>
#include <cmath>

namespace gapp::selection
{
    /* Calculate the cumulative distribution function of the population from the selection weights. */
    static std::vector<double> weightsToCdf(std::span<const double> weights)
    {
        GAPP_ASSERT(!weights.empty());
        GAPP_ASSERT(std::all_of(weights.begin(), weights.end(), detail::between(0.0, math::large<double>)));

        std::vector<double> cdf(weights.size());

        const double wmax = std::max(math::small<double>, detail::max_value(weights.begin(), weights.end()));
        const double wsum = std::transform_reduce(weights.begin(), weights.end(), 0.0, std::plus{}, detail::divide_by(wmax));
        const double idiv = std::isfinite(1.0 / wsum) ? (1.0 / wsum) : 1.0;
        const double corr = std::isfinite(1.0 / wsum) ? 0.0 : (1.0 / weights.size());

        std::transform_inclusive_scan(weights.begin(), weights.end(), cdf.begin(), std::plus{}, detail::multiply_add(idiv / wmax, corr));

        return cdf;
    }


    void Roulette::prepareSelectionsImpl(const GaInfo&, const PopulationView& pop)
    {
        GAPP_ASSERT(!pop.empty());

        FitnessVector fvec = detail::toFitnessVector(pop);

        /* Adjust for negative fitness values. */
        const double lowest = detail::min_value(fvec.begin(), fvec.end());
        const double offset = (-2.0 / 3.0) * lowest;
        const double scale  = 1.0 / 3.0;

        if (lowest < 0.0) std::transform(fvec.begin(), fvec.end(), fvec.begin(), detail::multiply_add(scale, offset));

        cdf_ = weightsToCdf(fvec);
    }

    size_t Roulette::selectImpl(const GaInfo&, const PopulationView&) const
    {
        return rng::sampleCdf(cdf_);
    }


    size_t Tournament::selectImpl(const GaInfo&, const PopulationView& pop) const
    {
        GAPP_ASSERT(!pop.empty());

        small_vector<size_t> candidate_indexes(tourney_size_);
        std::generate(candidate_indexes.begin(), candidate_indexes.end(), [&] { return rng::randomIndex(pop); });

        return *detail::max_element(candidate_indexes.begin(), candidate_indexes.end(), [&](size_t idx) { return pop[idx].fitness[0]; });
    }


    Rank::Rank(NonNegative<double> min_weight, NonNegative<double> max_weight) noexcept :
        min_weight_(min_weight), max_weight_(max_weight)
    {
        GAPP_ASSERT(min_weight <= max_weight, "The maximum selection weight can't be less than the minimum.");
    }

    void Rank::weights(NonNegative<double> min_weight, NonNegative<double> max_weight) noexcept
    {
        GAPP_ASSERT(min_weight <= max_weight, "The maximum selection weight can't be less than the minimum.");

        min_weight_ = min_weight;
        max_weight_ = max_weight;
    }

    void Rank::prepareSelectionsImpl(const GaInfo&, const PopulationView& pop)
    {
        GAPP_ASSERT(pop.size() > 1);

        FitnessVector fvec = detail::toFitnessVector(pop);
        const auto indices = detail::argsort(fvec.begin(), fvec.end());

        const double weight_step = (max_weight_ - min_weight_) / (pop.size() - 1);
        for (double weight = min_weight_; size_t idx : indices)
        {
            fvec[idx] = weight;
            weight += weight_step;
        }

        cdf_ = weightsToCdf(fvec);
    }

    size_t Rank::selectImpl(const GaInfo&, const PopulationView&) const
    {
        return rng::sampleCdf(cdf_);
    }


    void Sigma::prepareSelectionsImpl(const GaInfo&, const PopulationView& pop)
    {
        FitnessVector fvec = detail::toFitnessVector(pop);

        const double fmean = math::mean(fvec);
        const double fdev = math::stdDev(fvec, fmean);
        const double idiv = 1.0 / std::clamp(fdev, math::small<double>, math::large<double>);

        for (double& f : fvec)
        {
            const double weight = 1.0 + (f - fmean) * idiv;
            f = std::clamp(weight, 0.0, math::large<double>);
        }

        cdf_ = weightsToCdf(fvec);
    }

    size_t Sigma::selectImpl(const GaInfo&, const PopulationView&) const
    {
        return rng::sampleCdf(cdf_);
    }


    Boltzmann::Boltzmann(TemperatureFunction f) noexcept
    {
        GAPP_ASSERT(f, "The temperature function can't be a nullptr.");

        temperature_ = std::move(f);
    }

    void Boltzmann::prepareSelectionsImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(!pop.empty());
        GAPP_ASSERT(temperature_);

        FitnessVector fvec = detail::toFitnessVector(pop);

        const auto [fmin, fmax] = detail::minmax_value(fvec.begin(), fvec.end());
        const double df = std::clamp(fmax - fmin, math::small<double>, math::large<double>);
        const double temperature = std::invoke(temperature_, ga.generation_cntr(), ga.max_gen());

        for (double& f : fvec)
        {
            const double fnorm = (f / df) - (fmin / df);
            f = std::min(math::large<double>, std::exp(fnorm / temperature));
        }

        cdf_ = weightsToCdf(fvec);
    }

    size_t Boltzmann::selectImpl(const GaInfo&, const PopulationView&) const
    {
        return rng::sampleCdf(cdf_);
    }

    double Boltzmann::boltzmannDefaultTemp(size_t gen, size_t max_gen) noexcept
    {
        static constexpr double Tmin = 0.2;
        static constexpr double Tmax = 4.0;
        static constexpr double Tbeg = 3.0;
        static constexpr double Vtd = 10.0;

        return -Tmax / (1.0 + std::exp(-Vtd * (double(gen) / max_gen) + Tbeg)) + Tmax + Tmin;
    }


    Lambda::Lambda(SelectionCallable f) noexcept :
        selection_(std::move(f))
    {
        GAPP_ASSERT(selection_, "The selection method can't be a nullptr.");
    }

    size_t Lambda::selectImpl(const GaInfo& ga, const PopulationView& pop) const
    {
        GAPP_ASSERT(selection_);
        return selection_(ga, pop);
    }

} // namespace gapp::selection
