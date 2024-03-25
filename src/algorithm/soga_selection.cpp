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

        const double wmax = std::max(math::small<double>, *std::max_element(weights.begin(), weights.end()));
        const double wsum = std::transform_reduce(weights.begin(), weights.end(), 0.0, std::plus{}, detail::divide_by(wmax));
        const double idiv = std::isfinite(1.0 / wsum) ? (1.0 / wsum) : 1.0;
        const double corr = std::isfinite(1.0 / wsum) ? 0.0 : (1.0 / weights.size());

        std::transform_inclusive_scan(weights.begin(), weights.end(), cdf.begin(), std::plus{}, detail::multiply_add(idiv / wmax, corr));

        return cdf;
    }


    void Roulette::prepareSelectionsImpl(const GaInfo&, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(!fmat.empty());

        FitnessVector fvec = detail::toFitnessVector(fmat.begin(), fmat.end());

        /* Roulette selection wouldn't work for negative fitness values. */
        const double lowest = *std::min_element(fvec.begin(), fvec.end());
        const double offset = -(2.0 / 3.0) * std::min(0.0, lowest); // only adjust if there are negative fitness values
        const double scale = (lowest >= 0.0) ? 1.0 : (1.0 / 3.0);   // scaling to prevent overflow

        std::transform(fvec.begin(), fvec.end(), fvec.begin(), detail::multiply_add(scale, offset));

        cdf_ = weightsToCdf(fvec);
    }

    size_t Roulette::selectImpl(const GaInfo&, const FitnessMatrix&) const
    {
        return rng::sampleCdf<double>(cdf_);
    }


    size_t Tournament::selectImpl(const GaInfo&, const FitnessMatrix& fmat) const
    {
        GAPP_ASSERT(fmat.size() >= tourney_size_);
        GAPP_ASSERT(fmat.ncols() == 1);

        small_vector<size_t> candidate_indices(tourney_size_);
        std::generate(candidate_indices.begin(), candidate_indices.end(), [&]{ return rng::randomIndex(fmat); });

        return *std::max_element(candidate_indices.begin(), candidate_indices.end(),
        [&](size_t lhs, size_t rhs) noexcept
        {
            return fmat[lhs][0] < fmat[rhs][0];
        });
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

    void Rank::prepareSelectionsImpl(const GaInfo&, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(fmat.ncols() == 1);

        const auto indices = detail::argsort(fmat.begin(), fmat.end(), [](const auto& lhs, const auto& rhs) { return lhs[0] < rhs[0]; });

        std::vector<double> weights(fmat.size());
        for (size_t i = 0; i < indices.size(); i++)
        {
            const double t = i / (weights.size() - 1.0);
            weights[indices[i]] = std::lerp(min_weight_, max_weight_, t);
        }

        cdf_ = weightsToCdf(weights);
    }

    size_t Rank::selectImpl(const GaInfo&, const FitnessMatrix&) const
    {
        return rng::sampleCdf<double>(cdf_);
    }


    void Sigma::prepareSelectionsImpl(const GaInfo&, const FitnessMatrix& fmat)
    {
        FitnessVector fvec = detail::toFitnessVector(fmat.begin(), fmat.end());
        const double fmean = math::mean(fvec);
        const double fdev = math::stdDev(fvec, fmean);
        const double idiv = 1.0 / std::clamp(fdev, math::small<double>, math::large<double>);

        std::transform(fvec.begin(), fvec.end(), fvec.begin(), [&](double f)
        {
            const double weight = 1.0 + (f - fmean) * idiv;

            return std::clamp(weight, 0.0, math::large<double>);  /* If ( fitness < (f_mean - scale * SD) ) the weight could be negative. */
        });

        cdf_ = weightsToCdf(fvec);
    }

    size_t Sigma::selectImpl(const GaInfo&, const FitnessMatrix&) const
    {
        return rng::sampleCdf<double>(cdf_);
    }


    Boltzmann::Boltzmann(TemperatureFunction f) noexcept
    {
        GAPP_ASSERT(f, "The temperature function can't be a nullptr.");

        temperature_ = std::move(f);
    }

    void Boltzmann::prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(fmat.ncols() == 1);

        FitnessVector fvec = detail::toFitnessVector(fmat.begin(), fmat.end());
        const auto [fmin, fmax] = std::minmax_element(fvec.begin(), fvec.end());
        const double df = std::clamp(*fmax - *fmin, math::small<double>, math::large<double>);
        const double temperature = temperature_(ga.generation_cntr(), ga.max_gen());

        // can't capture the iterators by ref or value here
        std::transform(fvec.begin(), fvec.end(), fvec.begin(), [&, fmin = *fmin](double f) noexcept
        {
            const double fnorm = f / df - fmin / df; // normalize the fitness values to prevent overflows with std::exp

            return std::min(math::large<double>, std::exp(fnorm / temperature));
        });

        cdf_ = weightsToCdf(fvec);
    }

    size_t Boltzmann::selectImpl(const GaInfo&, const FitnessMatrix&) const
    {
        return rng::sampleCdf<double>(cdf_);
    }

    double Boltzmann::boltzmannDefaultTemp(size_t gen, size_t max_gen) noexcept
    {
        constexpr double Tmin = 0.2;
        constexpr double Tmax = 4.0;
        constexpr double tbeg = 3.0;
        constexpr double Vtd = 10.0;

        return -Tmax / (1.0 + std::exp(-Vtd * (double(gen) / max_gen) + tbeg)) + Tmax + Tmin;
    }


    Lambda::Lambda(SelectionCallable f) noexcept
    {
        GAPP_ASSERT(f, "The selection method can't be a nullptr.");

        selection_ = std::move(f);
    }

    size_t Lambda::selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const
    {
        GAPP_ASSERT(selection_);

        return selection_(ga, fmat);
    }

} // namespace gapp::selection