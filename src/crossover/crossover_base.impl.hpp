/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BASE_IMPL_HPP
#define GA_CROSSOVER_BASE_IMPL_HPP

#include "crossover_base.decl.hpp"
#include "../core/ga_base.hpp"
#include "../utility/rng.hpp"
#include <algorithm>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover
{
    template<Gene T>
    Crossover<T>::Crossover(double pc)
    {
        crossover_rate(pc);
    }

    template<Gene T>
    void Crossover<T>::crossover_rate(double pc)
    {
        if (!(0.0 <= pc && pc <= 1.0))
        {
            throw std::invalid_argument("The crossover probability must be in the closed range [0.0, 1.0]");
        }

        pc_ = pc;
    }

    template<Gene T>
    CandidatePair<T> Crossover<T>::operator()(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        assert(0.0 <= pc_ && pc_ <= 1.0);

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