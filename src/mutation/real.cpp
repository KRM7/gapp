/*
*  MIT License
*
*  Copyright (c) 2021 Krisztián Rugási
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in all
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

#include "real.hpp"
#include "../utils.h"
#include "../rng.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace genetic_algorithm::mutation::real
{
    void Uniform::mutate(const GA<double>& ga, Candidate<double>& candidate) const
    {
        GA_UNUSED(ga);

        if (candidate.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Uniform mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                candidate.chromosome[i] = rng::randomReal(bounds_[i].first, bounds_[i].second);
            }
        }
    }

    NonUniform::NonUniform(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double beta = 2.0) :
        BoundedMutation(bounds, pm)
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

    void NonUniform::mutate(const GA<double>& ga, Candidate<double>& candidate) const
    {
        if (candidate.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Non-Uniform mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            /* Perform mutation on the gene with pm probability. */
            if (rng::randomReal() < pm_)
            {
                double interval = bounds_[i].second - bounds_[i].first;
                double rand = rng::randomReal();
                double sign = rng::randomBool() ? 1.0 : -1.0;
                double exponent = std::pow(1.0 - double(ga.generation_cntr()) / ga.max_gen(), beta_);

                candidate.chromosome[i] += sign * interval * (1.0 - std::pow(rand, exponent));

                /* The mutated gene might be outside the allowed range. */
                candidate.chromosome[i] = std::clamp(candidate.chromosome[i], bounds_[i].first, bounds_[i].second);
            }
        }
    }

    Gauss::Gauss(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double sigma = 6.0) :
        BoundedMutation(bounds, pm)
    {
        this->sigma(sigma);
    }

    void Gauss::sigma(double sigma)
    {
        if (!(0.0 < sigma && sigma <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Sigma must be a nonnegative, finite value.");
        }

        sigma_ = sigma;
    }

    void Gauss::mutate(const GA<double>& ga, Candidate<double>& candidate) const
    {

    }

    Polynomial::Polynomial(const std::vector<std::pair<double, double>>& bounds, double pm = 0.01, double param = 6.0)
        : BoundedMutation(bounds, pm)
    {
        this->param(param);
    }

    void Polynomial::param(double param)
    {
        if (!(0.0 <= param && param <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Eta must be a nonnegative, finite value.");
        }

        param_ = param;
    }

    void Polynomial::mutate(const GA<double>& ga, Candidate<double>& candidate) const
    {
        GA_UNUSED(ga);

        if (candidate.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Polynomial mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                double rand = rng::randomReal();
                if (rand < 0.5)
                {
                    double delta = std::pow(2.0 * rand, 1.0 / (1.0 + param_)) - 1.0;
                    candidate.chromosome[i] += delta * (candidate.chromosome[i] - bounds_[i].first);
                }
                else
                {
                    double delta = 1.0 - std::pow(2.0 - 2.0 * rand, 1.0 / (1.0 + param_));
                    candidate.chromosome[i] += delta * (bounds_[i].second - candidate.chromosome[i]);
                }
                /* The mutated gene will always be in the allowed range. */
            }
        }
    }

    void Boundary::mutate(const GA<double>& ga, Candidate<double>& candidate) const
    {
        GA_UNUSED(ga);

        if (candidate.chromosome.size() != bounds_.size())
        {
            throw std::invalid_argument("The length of the chromosome must be the same as the bounds vector to perform the Boundary mutation.");
        }

        for (size_t i = 0; i < candidate.chromosome.size(); i++)
        {
            if (rng::randomReal() < pm_)
            {
                candidate.chromosome[i] = rng::randomBool() ? bounds_[i].first : bounds_[i].second;
            }
        }
    }

} // namespace genetic_algorithm::mutation::real