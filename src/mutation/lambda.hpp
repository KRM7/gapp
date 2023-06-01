/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_LAMBDA_HPP
#define GA_MUTATION_LAMBDA_HPP

#include "mutation_base.hpp"
#include "../population/candidate.hpp"
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
        using MutationCallable = std::function<void(const GA<T>&, const Candidate<T>&, Chromosome<T>&)>;

        constexpr explicit Lambda(MutationCallable f) noexcept :
            Mutation<T>(0.01)
        {
            GA_ASSERT(f, "The mutation method can't be a nullptr.");

            mutate_ = std::move(f);
        }

    private:
        MutationCallable mutate_;

        void mutate(const GA<T>& ga, const Candidate<T>& candidate, Chromosome<T>& chromosome) const override
        {
            GA_ASSERT(mutate_);

            mutate_(ga, candidate, chromosome);
        }
    };

} // namespace gapp::mutation

#endif // !GA_MUTATION_LAMBDA_HPP