/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SOGA_REPLACEMENT_HPP
#define GA_ALGORITHM_SOGA_REPLACEMENT_HPP

#include "replacement_base.hpp"
#include "../population/population.hpp"
#include <vector>
#include <functional>
#include <cstddef>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

namespace genetic_algorithm::replacement
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /**
    * A population update method that selects only the child Candidates from the
    * combined parent and child populations, and uses these as the population of the next generation of the algorithm. \n
    * 
    * If the number of children is greater than the population size used in the algorithm,
    * only the first pop_size children will be selected.
    */
    class KeepChildren final : public Replacement
    {
    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last) override;
    };


    /**
    * A population update method that selects the candidates of the next generation using elitism. \n
    * 
    * Of the combined parent and child populations, the N best candidates of the parent
    * population are carried over to the next population, while the remaining
    * (pop_size - N) slots are filled by the first (pop_size - N) child solutions. \n
    * 
    * If N = 0, this is equivalent to only keeping the children for the next generation (KeepChildren).
    */
    class Elitism final : public Replacement
    {
    public:
        /**
        * Create an elitist population update operator.
        * 
        * @param n The number of solutions from the population that will be copied to the next generation of the algorithm.
        */
        constexpr Elitism(size_t n = 1) noexcept :
            n_(n)
        {}

        /**
        * Set the number of elite solutions used to n.
        * 
        * @param n The number of solutions from the population that will be copied to the next generation of the algorithm.
        */
        constexpr void elite_num(size_t n) noexcept { n_ = n; }

        /** @returns The number of elite solutions used. */
        [[nodiscard]]
        constexpr size_t elite_num() const noexcept { return n_; }

    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last) override;

        size_t n_;
    };


    /**
    * A population update method that selects the best (pop_size) candidates of the combined
    * parent and child populations, and uses these as the candidates of the next generation of the algorithm. \n
    */
    class KeepBest final : public Replacement
    {
    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last) override;
    };


    /*
    * Wraps a callable with the right signature so that it can be used as a population replacement
    * method in the single-objective GAs.
    */
    class Lambda final : public Replacement
    {
    public:
        using ReplacementCallable = std::function<std::vector<size_t>(const GaInfo&, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator)>;

        explicit Lambda(ReplacementCallable f) noexcept;

    private:
        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last) override;

        ReplacementCallable replacement_;
    };

} // namespace genetic_algorithm::replacement


#endif // !GA_ALGORITHM_SOGA_REPLACEMENT_HPP