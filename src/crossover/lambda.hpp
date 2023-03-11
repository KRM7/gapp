/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_LAMBDA_HPP
#define GA_CROSSOVER_LAMBDA_HPP

#include "crossover_base.hpp"
#include "../population/candidate.hpp"
#include "../utility/utility.hpp"
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
        using CrossoverCallable = std::function<CandidatePair<T>(const GA<T>&, const Candidate<T>&, const Candidate<T>&)>;

        constexpr explicit Lambda(CrossoverCallable f) noexcept
        {
            GA_ASSERT(f, "The crossover method can't be a nullptr.");

            crossover_ = std::move(f);
        }

    private:
        CrossoverCallable crossover_;

        CandidatePair<T> crossover(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const override
        {
            GA_ASSERT(crossover_);

            return crossover_(ga, parent1, parent2);
        }
    };

} // namespace genetic_algorithm::crossover

#endif