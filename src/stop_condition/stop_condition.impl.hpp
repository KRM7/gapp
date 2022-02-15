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

#ifndef GA_STOP_CONDITION_IMPL_HPP
#define GA_STOP_CONDITION_IMPL_HPP

#include "stop_condition.hpp"
#include "../base_ga.h"
#include "../population.hpp"
#include "../math.hpp"
#include "../utils.hpp"

#include <algorithm>
#include <stdexcept>

namespace genetic_algorithm::stopping
{
    template<gene T>
    FitnessEvals<T>::FitnessEvals(size_t max_fitness_evals)
        : StopCondition<T>()
    {
        this->max_fitness_evals(max_fitness_evals);
    }

    template<gene T>
    FitnessEvals<T>::FitnessEvals(const GA<T>& ga, size_t max_fitness_evals)
        : StopCondition<T>()
    {
        GA_UNUSED(ga);

        this->max_fitness_evals(max_fitness_evals);
    }

    template<gene T>
    void FitnessEvals<T>::max_fitness_evals(size_t max_fitness_evals)
    {
        max_fitness_evals_ = max_fitness_evals;
    }

    template<gene T>
    bool FitnessEvals<T>::operator()(const GA<T>& ga)
    {
        return (ga.num_fitness_evals() >= max_fitness_evals_);
    }


    template<gene T>
    FitnessValue<T>::FitnessValue(const std::vector<double>& fitness_threshold) :
        StopCondition<T>()
    {
        this->fitness_threshold(fitness_threshold);
    }

    template<gene T>
    FitnessValue<T>::FitnessValue(const GA<T>& ga, const std::vector<double>& fitness_threshold) :
        StopCondition<T>()
    {
        GA_UNUSED(ga);

        this->fitness_threshold(fitness_threshold);
    }

    template<gene T>
    void FitnessValue<T>::fitness_threshold(const std::vector<double>& fitness_threshold)
    {
        if (fitness_threshold.empty())
        {
            throw std::invalid_argument("Empty fitness threshold vector.");
        }

        fitness_threshold_ = fitness_threshold;
    }

    template<gene T>
    bool FitnessValue<T>::operator()(const GA<T>& ga)
    {
        if (ga.num_objectives() != fitness_threshold_.size())
        {
            throw std::domain_error("The size of the fitness threshold vector does not match the size of the fitness vectors.");
        }

        return std::any_of(ga.population().begin(), ga.population().end(),
        [this](const Candidate<T>& sol)
        {
            return detail::paretoCompareLess(fitness_threshold_, sol.fitness);
        });
    }


    template<gene T>
    FitnessMeanStall<T>::FitnessMeanStall(size_t patience, double delta)
        : StopCondition<T>()
    {
        this->patience(patience);
        this->delta(delta);
    }

    template<gene T>
    FitnessMeanStall<T>::FitnessMeanStall(const GA<T>& ga, size_t patience, double delta)
        : StopCondition<T>()
    {
        GA_UNUSED(ga);

        this->patience(patience);
        this->delta(delta);
    }

    template<gene T>
    void FitnessMeanStall<T>::patience(size_t patience)
    {
        patience_ = patience;
        resetCntr();
    }

    template<gene T>
    void FitnessMeanStall<T>::delta(double delta)
    {
        delta_ = delta;
    }

    template<gene T>
    void FitnessMeanStall<T>::resetCntr()
    {
        cntr_ = patience_ + 1;
    }

    template<gene T>
    bool FitnessMeanStall<T>::operator()(const GA<T>& ga)
    {
        auto current_mean = populationFitnessMean(ga.population());

        /* Init on first gen. */
        if (ga.generation_cntr() == 0)
        {
            resetCntr();
            best_fitness_mean_ = current_mean;

            return false;
        }

        bool improved = false;
        for (size_t i = 0; i < current_mean.size(); i++)
        {
            if (current_mean[i] >= best_fitness_mean_[i] + delta_)
            {
                best_fitness_mean_[i] = current_mean[i];
                improved = true;
            }
        }

        if (improved)
        {
            resetCntr();
        }
        else
        {
            --cntr_;
        }

        return cntr_ == 0;
    }


    template<gene T>
    FitnessBestStall<T>::FitnessBestStall(size_t patience, double delta) :
        StopCondition<T>()
    {
        this->patience(patience);
        this->delta(delta);
    }

    template<gene T>
    FitnessBestStall<T>::FitnessBestStall(const GA<T>& ga, size_t patience, double delta) :
        StopCondition<T>()
    {
        GA_UNUSED(ga);

        this->patience(patience);
        this->delta(delta);
    }

    template<gene T>
    void FitnessBestStall<T>::patience(size_t patience)
    {
        patience_ = patience;
        resetCntr();
    }

    template<gene T>
    void FitnessBestStall<T>::delta(double delta)
    {
        delta_ = delta;
    }

    template<gene T>
    void FitnessBestStall<T>::resetCntr()
    {
        cntr_ = patience_ + 1;
    }

    template<gene T>
    bool FitnessBestStall<T>::operator()(const GA<T>& ga)
    {
        auto current_max = populationFitnessMax(ga.population());

        /* Init on first gen. */
        if (ga.generation_cntr() == 0)
        {
            resetCntr();
            best_fitness_max_ = current_max;

            return false;
        }

        bool improved = false;
        for (size_t i = 0; i < current_max.size(); i++)
        {
            if (current_max[i] >= best_fitness_max_[i] + delta_)
            {
                best_fitness_max_[i] = current_max[i];
                improved = true;
            }
        }

        if (improved)
        {
            resetCntr();
        }
        else
        {
            --cntr_;
        }

        return cntr_ == 0;
    }

} // namespace genetic_algorithm::stopping

#endif // !GA_STOP_CONDITION_IMPL_HPP