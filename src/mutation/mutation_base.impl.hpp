/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_IMPL_HPP
#define GA_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../core/ga_base.decl.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <utility>
#include <stdexcept>

namespace genetic_algorithm::mutation
{
    template<Gene T>
    void Mutation<T>::operator()(const GA<T>& ga, Candidate<T>& candidate) const
    {
        GA_ASSERT(!candidate.is_evaluated || candidate.fitness.size() == ga.num_objectives());
        GA_ASSERT(ga.variable_chrom_len() || candidate.chromosome.size() == ga.chrom_len());

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
            thread_local Chromosome<T> old_chromosome; old_chromosome = candidate.chromosome;

            /* Assume that the call to mutate() doesn't change the fitness vector. */
            mutate(ga, candidate);
            candidate.is_evaluated = (candidate.chromosome == old_chromosome);
        }

        if (!ga.variable_chrom_len() && candidate.chromosome.size() != ga.chrom_len())
        {
            GA_THROW(std::logic_error, "The mutation resulted in a candidate with incorrect chromosome length.");
        }
    }

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_IMPL_HPP