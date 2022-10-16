/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_LAMBDA_HPP
#define GA_CROSSOVER_LAMBDA_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include <functional>
#include <utility>

namespace genetic_algorithm::crossover
{
    /*
    * Wraps a callable with the right signature so that it can be used as a crossover
    * method in the GAs.
    */
    template<Gene T>
    class Lambda final : public Crossover<T>
    {
    public:
        using CrossoverFunction = std::function<CandidatePair<T>(const GA<T>&, const Candidate<T>&, const Candidate<T>&)>;

        explicit Lambda(CrossoverFunction f)
            : Crossover<T>(), crossover_(std::move(f))
        {}

    private:
        CrossoverFunction crossover_;

        CandidatePair<T> crossover(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const override
        {
            return crossover_(ga, parent1, parent2);
        }
    };

} // namespace genetic_algorithm::crossover

#endif