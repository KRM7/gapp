/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "real.hpp"
#include "../algorithms/ga_info.hpp"
#include "../algorithms/real_ga.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace genetic_algorithm::mutation::real
{
    void Uniform::mutate(const GaInfo& ga, Candidate<double>& candidate) const
    {
        auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Uniform mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                candidate.chromosome[i] = rng::randomReal(bounds[i].first, bounds[i].second);
            }
        }
    }

    NonUniform::NonUniform(double pm, double beta) :
        Mutation(pm)
    {
        this->beta(beta);
    }

    void NonUniform::beta(double beta)
    {
        if (!(0.0 <= beta && beta <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("The beta parameter must be a nonnegative, finite value.");
        }

        beta_ = beta;
    }

    void NonUniform::mutate(const GaInfo& ga, Candidate<double>& candidate) const
    {
        auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Non-Uniform mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                double rand = rng::randomReal();
                double exponent = std::pow(1.0 - double(ga.generation_cntr()) / ga.max_gen(), beta_);

                if (rng::randomBool())
                {
                    candidate.chromosome[i] += (bounds[i].second - candidate.chromosome[i]) * (1.0 - std::pow(rand, exponent));
                }
                else
                {
                    candidate.chromosome[i] -= (candidate.chromosome[i] - bounds[i].first) * (1.0 - std::pow(rand, exponent));
                }
                /* The mutated gene might be outside the allowed range. */
                candidate.chromosome[i] = std::clamp(candidate.chromosome[i], bounds[i].first, bounds[i].second);
            }
        }
    }

    Gauss::Gauss(double pm, double sigma) :
        Mutation(pm)
    {
        this->sigma(sigma);
    }

    void Gauss::sigma(double sigma)
    {
        if (!(0.0 < sigma && sigma <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Sigma must be a positive, finite value.");
        }

        sigma_ = sigma;
    }

    void Gauss::mutate(const GaInfo& ga, Candidate<double>& candidate) const
    {
        auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Gauss mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() <= pm_)
            {
                double SD = (bounds[i].second - bounds[i].first) / sigma_;

                double delta = 0.0;
                for (size_t j = 0; j < 5; j++)
                {
                    delta = rng::randomNormal(0.0, SD);
                    if (bounds[i].first <= (candidate.chromosome[i] + delta) && (candidate.chromosome[i] + delta) <= bounds[i].second)
                    {
                        break;
                    }
                }
                candidate.chromosome[i] = std::clamp(candidate.chromosome[i] + delta, bounds[i].first, bounds[i].second);
            }
        }
    }

    Polynomial::Polynomial(double pm, double eta)
        : Mutation(pm)
    {
        this->eta(eta);
    }

    void Polynomial::eta(double eta)
    {
        if (!(0.0 <= eta && eta <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        eta_ = eta;
    }

    void Polynomial::mutate(const GaInfo& ga, Candidate<double>& candidate) const
    {
        auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Polynomial mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                double alpha = rng::randomReal();
                if (alpha <= 0.5)
                {
                    double delta = std::pow(2.0 * alpha, -(1.0 + eta_)) - 1.0;
                    candidate.chromosome[i] += delta * (candidate.chromosome[i] - bounds[i].first);
                }
                else
                {
                    double delta = 1.0 - std::pow(2.0 - 2.0 * alpha, -(1.0 + eta_));
                    candidate.chromosome[i] += delta * (bounds[i].second - candidate.chromosome[i]);
                }
                /* The mutated gene might be outside the allowed range. */
                candidate.chromosome[i] = std::clamp(candidate.chromosome[i], bounds[i].first, bounds[i].second);
            }
        }
    }

    void Boundary::mutate(const GaInfo& ga, Candidate<double>& candidate) const
    {
        auto& bounds = dynamic_cast<const RCGA&>(ga).limits();

        if (candidate.chromosome.size() != bounds.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Boundary mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                candidate.chromosome[i] = rng::randomBool() ? bounds[i].first : bounds[i].second;
            }
        }
    }

} // namespace genetic_algorithm::mutation::real