/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_IMPL_HPP
#define GA_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../core/ga_base.decl.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <utility>

namespace gapp::mutation
{
    template<typename T>
    void Mutation<T>::operator()(const GA<T>& ga, Candidate<T>& candidate) const
    {
        GAPP_ASSERT(candidate.fitness.empty() || candidate.fitness.size() == ga.num_objectives());
        GAPP_ASSERT(allow_variable_chrom_length() || candidate.chromosome.size() == ga.chrom_len());

        if (!candidate.is_evaluated())
        {
            mutate(ga, candidate, candidate.chromosome);
        }
        else
        {
            /* If the candidate is already evaluated (happens when the crossover didn't change the candidate),
               its current fitness vector is valid, and we can save a fitness function call when the mutation
               doesn't change the chromosome. */
            thread_local Candidate<T> old_candidate; old_candidate = candidate;

            mutate(ga, candidate, candidate.chromosome);
            if (candidate != old_candidate) candidate.fitness.clear();
        }

        GAPP_ASSERT(allow_variable_chrom_length() || candidate.chromosome.size() == ga.chrom_len(),
                  "The mutation resulted in a candidate with incorrect chromosome length.");
    }

} // namespace gapp::mutation

#endif // !GA_MUTATION_BASE_IMPL_HPP