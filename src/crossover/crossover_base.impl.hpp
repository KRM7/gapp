/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BASE_IMPL_HPP
#define GA_CROSSOVER_BASE_IMPL_HPP

#include "crossover_base.decl.hpp"
#include "../core/ga_base.decl.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <utility>

namespace gapp::crossover
{
    template<typename T>
    CandidatePair<T> Crossover<T>::operator()(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        GAPP_ASSERT(parent1.is_evaluated() && parent2.is_evaluated());
        GAPP_ASSERT(parent1.fitness.size() == parent2.fitness.size());
        GAPP_ASSERT(parent1.fitness.size() == ga.num_objectives());
        GAPP_ASSERT(allow_variable_chrom_length() || (parent1.chromosome.size() == ga.chrom_len() &&
                                                      parent2.chromosome.size() == ga.chrom_len()));

        /* Only need to perform the crossover with the set pc probability. Return early with (1 - pc) probability. */
        if (rng::randomReal() >= pc_)
        {
            return { parent1, parent2 };
        }

        /*
        * If the parents are the same, the crossover doesn't need to be performed.
        * This assumes that with 2 parents that have the same chromosomes, the children would be the same
        * as the parents, which is true for every crossover operator implemented, but
        * could be an issue for user defined crossovers.
        */
        if (parent1 == parent2)
        {
            return { parent1, parent2 };
        }

        /* Perform the actual crossover. */
        auto [child1, child2] = crossover(ga, parent1, parent2);

        GAPP_ASSERT(allow_variable_chrom_length() || child1.chromosome.size() == ga.chrom_len(),
                  "The crossover returned a candidate with incorrect chromosome length.");
        GAPP_ASSERT(allow_variable_chrom_length() || child2.chromosome.size() == ga.chrom_len(),
                  "The crossover returned a candidate with incorrect chromosome length.");

        child1.fitness.clear();
        child2.fitness.clear();

        /*
        * Check if either of the children are the same as one of the parents.
        * (This can happen in edge cases even if the parents are different.)
        * If one of the children are the same as one of the parents, then the fitness function
        * evaluation for that child can be skipped (if the fitness function is the same) by assigning
        * it the same fitness as the parent.
        */
        if (child1 == parent1)
        {
            child1.fitness = parent1.fitness;
        }
        else if (child1 == parent2)
        {
            child1.fitness = parent2.fitness;
        }
        if (child2 == parent1)
        {
            child2.fitness = parent1.fitness;
        }
        else if (child2 == parent2)
        {
            child2.fitness = parent2.fitness;
        }

        return { std::move(child1), std::move(child2) };
    }

} // namespace gapp::crossover

#endif // !GA_CROSSOVER_BASE_IMPL_HPP