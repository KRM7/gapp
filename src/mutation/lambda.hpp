/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_LAMBDA_HPP
#define GA_MUTATION_LAMBDA_HPP

#include "mutation_base.decl.hpp"
#include "../population/candidate.hpp"
#include <functional>
#include <utility>

namespace genetic_algorithm::mutation
{
    /*
    * Wraps a callable with the right signature so that it can be used as a mutation
    * method in the GAs.
    */
    template<Gene T>
    class Lambda final : public Mutation<T>
    {
    public:
        using MutationFunction = std::function<void(const GA<T>&, Candidate<T>&)>;

        explicit Lambda(MutationFunction f)
            : Mutation<T>(0.01), mutate_(std::move(f))
        {}

    private:
        MutationFunction mutate_;

        void mutate(const GA<T>& ga, Candidate<T>& candidate) const override
        {
            mutate_(ga, candidate);
        }
    };

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_LAMBDA_HPP