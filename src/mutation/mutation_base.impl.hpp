/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_IMPL_HPP
#define GA_MUTATION_BASE_IMPL_HPP

#include "mutation_base.decl.hpp"
#include "../core/ga_base.hpp"
#include <algorithm>
#include <utility>
#include <cassert>
#include <stdexcept>

namespace genetic_algorithm::mutation
{
    template<Gene T>
    void Mutation<T>::operator()(const GaInfo& ga, Candidate<T>& candidate) const
    {
        assert(!candidate.is_evaluated || candidate.fitness.size() == ga.num_objectives());
        assert(ga.variable_chrom_len() || candidate.chromosome.size() == ga.chrom_len());

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
            const auto old_chromosome = candidate.chromosome;
            auto old_fitness = candidate.fitness;

            mutate(ga, candidate);
            candidate.is_evaluated = false;

            if (candidate.chromosome == old_chromosome)
            {
                candidate.fitness = std::move(old_fitness);
                candidate.is_evaluated = true;
            }
        }

        if (!ga.variable_chrom_len() && candidate.chromosome.size() != ga.chrom_len())
        {
            GA_THROW(std::logic_error, "The mutation resulted in a candidate with incorrect chromosome length.");
        }
    }

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_IMPL_HPP