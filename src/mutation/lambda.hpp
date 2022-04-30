/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_LAMBDA_HPP
#define GA_MUTATION_LAMBDA_HPP

#include "mutation_base.decl.hpp"
#include "../population/candidate.hpp"
#include <functional>

namespace genetic_algorithm::mutation::dtl
{
    template<Gene T>
    class Lambda final : public Mutation<T>
    {
    public:
        using MutationFunction = std::function<void(const GaInfo&, Candidate<T>&)>;

        explicit Lambda(MutationFunction f);

    private:
        void mutate(const GaInfo& ga, Candidate<T>& candidate) const override;

        MutationFunction mutate_;
    };

} // namespace genetic_algorithm::mutation::dtl


/* IMPLEMENTATION */

#include <utility>
#include <cassert>

namespace genetic_algorithm::mutation::dtl
{
    template<Gene T>
    Lambda<T>::Lambda(MutationFunction f)
        : Mutation<T>(0.01), mutate_(std::move(f))
    {
    }

    template<Gene T>
    void Lambda<T>::mutate(const GaInfo& ga, Candidate<T>& candidate) const
    {
        assert(mutate_);
        mutate_(ga, candidate);
    }

} // namespace genetic_algorithm::mutation::dtl

#endif // !GA_MUTATION_LAMBDA_HPP