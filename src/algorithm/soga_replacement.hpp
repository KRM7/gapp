/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_REPLACEMENT_HPP
#define GA_ALGORITHM_SOGA_REPLACEMENT_HPP

#include "replacement_base.hpp"
#include "../core/population.hpp"
#include <vector>
#include <functional>
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

namespace gapp::replacement
{
    /**
    * A population update method that selects only the child candidates from the
    * combined parent and child populations, and uses these as the population of the
    * next generation of the %GA.
    * 
    * If the number of children is greater than the population size used in the algorithm,
    * only the first @p population_size children will be selected.
    */
    class KeepChildren final : public Replacement
    {
    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
    };


    /**
    * A population update method that selects the candidates of the next generation using elitism.
    * 
    * The operator has a single parameter @p N, which determines the number of candidates
    * that will be selected from the parent population. Of the combined parent and child populations,
    * the N best candidates of the parent population will be copied over to the next population,
    * while the remaining (pop_size - N) slots are filled by the first (pop_size - N) child solutions.
    * 
    * If N is equal to 0, this is equivalent to only keeping the children for the next generation
    * (ie. KeepChildren).
    */
    class Elitism final : public Replacement
    {
    public:
        /**
        * Create an elitist population update operator.
        * 
        * @param n The number of solutions from the parent population that will be
        *   copied to the next generation of the algorithm.
        */
        constexpr Elitism(size_t n = 1) noexcept :
            n_(n)
        {}

        /**
        * Set the number of elite solutions used.
        * 
        * @param n The number of solutions from the parent population that will be
        *   copied to the next generation of the algorithm.
        */
        constexpr void elite_num(size_t n) noexcept { n_ = n; }

        /** @returns The number of elite solutions used. */
        [[nodiscard]]
        constexpr size_t elite_num() const noexcept { return n_; }

    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;

        size_t n_;
    };


    /**
    * A population update method that selects the best @p pop_size candidates of the combined
    * parent and child populations, and uses these as the candidates of the next generation of the
    * algorithm.
    * 
    * The operator assumes fitness maximization.
    */
    class KeepBest final : public Replacement
    {
    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
    };


    /*
    * Wraps a callable with the right signature so that it can be used as a population replacement
    * policy in the single-objective algorithm.
    */
    class Lambda final : public Replacement
    {
    public:
        using ReplacementCallable = std::function<std::vector<size_t>(const GaInfo&, const FitnessMatrix&)>;

        explicit Lambda(ReplacementCallable f) noexcept;

    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;

        ReplacementCallable replacement_;
    };

} // namespace gapp::replacement


#endif // !GA_ALGORITHM_SOGA_REPLACEMENT_HPP