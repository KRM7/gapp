/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BASE_IMPL_HPP
#define GA_CROSSOVER_BASE_IMPL_HPP

#include "crossover_base.decl.hpp"
#include "../core/ga_base.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <utility>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover
{
    template<Gene T>
    CandidatePair<T> Crossover<T>::operator()(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        assert(0.0 <= pc_ && pc_ <= 1.0);
        assert(parent1.is_evaluated && parent2.is_evaluated);
        assert(parent1.fitness.size() == ga.num_objectives() && parent2.fitness.size() == ga.num_objectives());
        assert(ga.variable_chromosome_length() || (parent1.chromosome.size() == ga.chrom_len() &&
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

        if (!ga.variable_chromosome_length())
        {
            if (child1.chromosome.size() != ga.chrom_len() || child2.chromosome.size() != ga.chrom_len())
            {
                GA_THROW(std::logic_error, "The crossover function returned a candidate with incorrect chromosome length.");
            }
        }

        child1.is_evaluated = false;
        child2.is_evaluated = false;

        /*
        * Check if either of the children are the same as one of the parents.
        * (This can happen in edge cases even if the parents are different.)
        * If one of the children are the same as one of the parents, then the fitness function
        * evaluation for that child can be skipped (if the fitness function is the same) by assigning
        * it the same fitness as the parent.
        */
        if (child1.chromosome == parent1.chromosome)
        {
            child1.fitness = parent1.fitness;
            child1.is_evaluated = parent1.is_evaluated;
        }
        else if (child1.chromosome == parent2.chromosome)
        {
            child1.fitness = parent2.fitness;
            child1.is_evaluated = parent2.is_evaluated;
        }
        if (child2.chromosome == parent1.chromosome)
        {
            child2.fitness = parent1.fitness;
            child2.is_evaluated = parent1.is_evaluated;
        }
        else if (child2.chromosome == parent2.chromosome)
        {
            child2.fitness = parent2.fitness;
            child2.is_evaluated = parent2.is_evaluated;
        }

        return { std::move(child1), std::move(child2) };
    }

} // namespace genetic_algorithm::crossover

#endif // !GA_CROSSOVER_BASE_IMPL_HPP