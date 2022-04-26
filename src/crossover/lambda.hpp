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
        using CrossoverFunction_t = std::function<CandidatePair<T>(const GaInfo&, const Candidate<T>&, const Candidate<T>&)>;

        Lambda(double pc, CrossoverFunction_t f);
        Lambda(double pc, std::function<CandidatePair<T>(const Candidate<T>&, const Candidate<T>&)> f);

    private:
        CandidatePair<T> crossover(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const override;

        CrossoverFunction_t crossover_;
    };

} // namespace genetic_algorithm::crossover::dtl


#include <utility>
#include <stdexcept>
#include <cassert>

namespace genetic_algorithm::crossover::dtl
{
    template<Gene T>
    Lambda<T>::Lambda(double pc, CrossoverFunction_t f)
    {
        if (!f)
        {
            throw std::invalid_argument("The crossover function can't be a nullptr.");
        }

        crossover_ = std::move(f);
    }

    template<Gene T>
    Lambda<T>::Lambda(double pc, std::function<CandidatePair<T>(const Candidate<T>&, const Candidate<T>&)> f)
    {
        if (!f)
        {
            throw std::invalid_argument("The crossover function can't be a nullptr.");
        }

        crossover_ = [f = std::move(f)](const GaInfo&, const Candidate<T>& parent1, const Candidate<T>& parent2)
        {
            return f(parent1, parent2);
        }
    }

    template<Gene T>
    CandidatePair<T> Lambda<T>::crossover(const GaInfo& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const
    {
        assert(crossover_);

        return crossover_(ga, parent1, parent2);
    }

} // namespace genetic_algorithm::crossover::dtl

#endif