/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_LAMBDA_HPP
#define GA_CROSSOVER_LAMBDA_HPP

#include "crossover_base.decl.hpp"
#include "../population/candidate.hpp"
#include <functional>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    class Lambda final : public Crossover<T>
    {
    public:
        using CrossoverFunction = std::function<CandidatePair<T>(const GaInfo&, const Candidate<T>&, const Candidate<T>&)>;

        explicit Lambda(CrossoverFunction f);

    private:
        CandidatePair<T> crossover(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const override;

        CrossoverFunction crossover_;
    };

} // namespace genetic_algorithm::crossover::dtl


/* IMPLEMENTATION */

#include <utility>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    Lambda<T>::Lambda(CrossoverFunction f)
        : Crossover<T>(), crossover_(std::move(f))
    {}

    template<Gene T>
    CandidatePair<T> Lambda<T>::crossover(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        return crossover_(ga, parent1, parent2);
    }

} // namespace genetic_algorithm::crossover::dtl

#endif