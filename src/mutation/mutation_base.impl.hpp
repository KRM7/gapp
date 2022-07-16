/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_IMPL_HPP
#define GA_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../core/ga_base.hpp"
#include <algorithm>
#include <stdexcept>

namespace genetic_algorithm::mutation
{
    template<Gene T>
    Mutation<T>::Mutation(double pm)
    {
        mutation_rate(pm);
    }

    template<Gene T>
    void Mutation<T>::mutation_rate(double pm)
    {
        if (!(0.0 <= pm && pm <= 1.0))
        {
            throw std::invalid_argument("The mutation probability must be in the closed range [0.0, 1.0]");
        }

        pm_ = pm;
    }

    template<Gene T>
    void Mutation<T>::operator()(const GaInfo& ga, Candidate<T>& candidate) const
    {
        if (!candidate.is_evaluated)
        {
            mutate(ga, candidate);
            candidate.is_evaluated = false;
        }
        else
        {
            /* If the candidate is already evaluated (happens when the crossover didn't change the candidate),
               its current fitness vector is valid, and we can save a fitness function call when the mutation
               doesn't change the chromosome. */
            auto old_chromosome = candidate.chromosome;

            #ifndef NDEBUG
                auto old_fitness = candidate.fitness;
                mutate(ga, candidate);
                assert(candidate.fitness == old_fitness);
            #else
                mutate(ga, candidate);
            #endif

            /* Assume that the fitness vector of candidate was not changed by mutate() (it shouldn't have been). */
            candidate.is_evaluated = (candidate.chromosome == old_chromosome);
        }
    }

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_IMPL_HPP