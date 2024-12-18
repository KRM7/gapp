/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CROSSOVER_LAMBDA_HPP
#define GAPP_CROSSOVER_LAMBDA_HPP

#include "crossover_base.hpp"
#include "../core/candidate.hpp"
#include "../utility/utility.hpp"
#include <functional>
#include <utility>

namespace gapp::crossover
{
    /*
    * Wraps a callable with the right signature so that it can be used as a crossover
    * operator in the GAs.
    */
    template<typename T>
    class Lambda final : public Crossover<T>
    {
    public:
        using CrossoverCallable = std::function<CandidatePair<T>(const GaInfo&, const Candidate<T>&, const Candidate<T>&)>;

        constexpr explicit Lambda(CrossoverCallable f) noexcept :
            crossover_(std::move(f))
        {
            GAPP_ASSERT(crossover_, "The crossover method can't be a nullptr.");
        }

    private:
        CrossoverCallable crossover_;

        CandidatePair<T> crossover(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const override
        {
            GAPP_ASSERT(crossover_);
            return crossover_(ga, parent1, parent2);
        }
    };

} // namespace gapp::crossover

#endif // !GAPP_CROSSOVER_LAMBDA_HPP
