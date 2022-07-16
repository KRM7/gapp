/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../core/ga_info.hpp"
#include "../encoding/real.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace genetic_algorithm::mutation::real
{
    void Uniform::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        const auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Uniform mutation.");
        }

        size_t mutate_count = rng::randomBinomialApprox(candidate.chromosome.size(), mutation_rate());
        auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            candidate.chromosome[idx] = rng::randomReal(bounds[idx].first, bounds[idx].second);
        }
    }

    NonUniform::NonUniform(double pm, GeneType beta) :
        Mutation(pm)
    {
        this->beta(beta);
    }

    void NonUniform::beta(GeneType beta)
    {
        if (!(0.0 <= beta && beta <= std::numeric_limits<GeneType>::max()))
        {
            throw std::invalid_argument("The beta parameter must be a nonnegative, finite value.");
        }

        beta_ = beta;
    }

    void NonUniform::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        const auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Non-Uniform mutation.");
        }

        size_t mutate_count = rng::randomBinomialApprox(candidate.chromosome.size(), mutation_rate());
        auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            GeneType rand = rng::randomReal<GeneType>();
            GeneType exponent = std::pow(1.0 - GeneType(ga.generation_cntr()) / ga.max_gen(), beta_);

            GeneType multiplier = 1.0 - std::pow(rand, exponent);
            GeneType bound = rng::randomBool() ? bounds[idx].second : bounds[idx].first;

            candidate.chromosome[idx] += (bound - candidate.chromosome[idx]) * multiplier;
            candidate.chromosome[idx] = std::clamp(candidate.chromosome[idx], bounds[idx].first, bounds[idx].second);
        }
    }

    Gauss::Gauss(double pm, GeneType sigma) :
        Mutation(pm)
    {
        this->sigma(sigma);
    }

    void Gauss::sigma(GeneType sigma)
    {
        if (!(0.0 < sigma && sigma <= std::numeric_limits<GeneType>::max()))
        {
            throw std::invalid_argument("Sigma must be a positive, finite value.");
        }

        sigma_ = sigma;
    }

    void Gauss::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        const auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Gauss mutation.");
        }

        size_t mutate_count = rng::randomBinomialApprox(candidate.chromosome.size(), mutation_rate());
        auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            GeneType SD = (bounds[idx].second - bounds[idx].first) / sigma_;

            candidate.chromosome[idx] += rng::randomNormal<GeneType>(0.0, SD);
            candidate.chromosome[idx] = std::clamp(candidate.chromosome[idx], bounds[idx].first, bounds[idx].second);
        }
    }

    Polynomial::Polynomial(double pm, GeneType eta)
        : Mutation(pm)
    {
        this->eta(eta);
    }

    void Polynomial::eta(GeneType eta)
    {
        if (!(0.0 <= eta && eta <= std::numeric_limits<GeneType>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        eta_ = eta;
    }

    void Polynomial::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        const auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Polynomial mutation.");
        }

        size_t mutate_count = rng::randomBinomialApprox(candidate.chromosome.size(), mutation_rate());
        auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            GeneType alpha = rng::randomReal<GeneType>();
            if (alpha <= 0.5)
            {
                GeneType delta = std::pow(2.0 * alpha, -(1.0 + eta_)) - 1.0;
                candidate.chromosome[idx] += delta * (candidate.chromosome[idx] - bounds[idx].first);
            }
            else
            {
                GeneType delta = 1.0 - std::pow(2.0 - 2.0 * alpha, -(1.0 + eta_));
                candidate.chromosome[idx] += delta * (bounds[idx].second - candidate.chromosome[idx]);
            }
            candidate.chromosome[idx] = std::clamp(candidate.chromosome[idx], bounds[idx].first, bounds[idx].second);
        }
    }

    void Boundary::mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const
    {
        const auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Boundary mutation.");
        }

        size_t mutate_count = rng::randomBinomialApprox(candidate.chromosome.size(), mutation_rate());
        auto mutated_indices = rng::sampleUnique(0_sz, candidate.chromosome.size(), mutate_count);

        for (const auto& idx : mutated_indices)
        {
            candidate.chromosome[idx] = rng::randomBool() ? bounds[idx].first : bounds[idx].second;
        }
    }

} // namespace genetic_algorithm::mutation::real