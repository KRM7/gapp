/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_MUTATION_LAMBDA_HPP
#define GAPP_MUTATION_LAMBDA_HPP

#include "mutation_base.hpp"
#include "../core/candidate.hpp"
#include "../utility/utility.hpp"
#include <functional>
#include <utility>

namespace gapp::mutation
{
    /*
    * Wraps a callable with the right signature so that it can be used as a mutation
    * method in the GAs.
    */
    template<typename T>
    class Lambda final : public Mutation<T>
    {
    public:
        using MutationCallable = std::function<void(const GaInfo&, const Candidate<T>&, Chromosome<T>&)>;

        constexpr explicit Lambda(MutationCallable f) noexcept :
            Mutation<T>(0.01),
            mutate_(std::move(f))
        {
            GAPP_ASSERT(mutate_, "The mutation method can't be a nullptr.");
        }

    private:
        MutationCallable mutate_;

        void mutate(const GaInfo& ga, const Candidate<T>& candidate, Chromosome<T>& chromosome) const override
        {
            GAPP_ASSERT(mutate_);
            mutate_(ga, candidate, chromosome);
        }
    };

} // namespace gapp::mutation

#endif // !GAPP_MUTATION_LAMBDA_HPP
