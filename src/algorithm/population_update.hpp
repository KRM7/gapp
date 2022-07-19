/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_POPULATION_UPDATE_HPP
#define GA_ALGORITHM_POPULATION_UPDATE_HPP

#include "../population/population.hpp"
#include "../core/ga_info.hpp"
#include <concepts>
#include <vector>
#include <cstddef>

/** ... */
namespace genetic_algorithm::population_update
{
    using detail::FitnessVector;
    using detail::FitnessMatrix;

    /** 
     * Concept specifying the interface required for the population update methods.
     */
    template<typename T>
    concept Updater =
    requires(T updater, GaInfo ga, FitnessMatrix::const_iterator it)
    {
        requires std::copyable<T>;
        requires std::destructible<T>;

        { updater(ga, it, it, it) } -> std::same_as<std::vector<size_t>>;
    };

    /** ... */
    class KeepChildren
    {
    public:
        /** ... */
        std::vector<size_t> operator()(const GaInfo& ga,
                                       FitnessMatrix::const_iterator parents_first,
                                       FitnessMatrix::const_iterator children_first,
                                       FitnessMatrix::const_iterator children_last);
    };

    /** ... */
    class Elitism
    {
    public:
        /** ... */
        Elitism(size_t n = 1);

        /** ... */
        void elite_num(size_t n);

        /** ... */
        [[nodiscard]]
        size_t elite_num() noexcept { return n_; }

        /** ... */
        std::vector<size_t> operator()(const GaInfo& ga,
                                       FitnessMatrix::const_iterator parents_first,
                                       FitnessMatrix::const_iterator children_first,
                                       FitnessMatrix::const_iterator children_last);

    private:
        size_t n_;
    };

    /** ... */
    struct KeepBest
    {
    public:
        /** ... */
        std::vector<size_t> operator()(const GaInfo& ga,
                                       FitnessMatrix::const_iterator parents_first,
                                       FitnessMatrix::const_iterator children_first,
                                       FitnessMatrix::const_iterator children_last);
    };

    static_assert(Updater<KeepChildren>);
    static_assert(Updater<Elitism>);
    static_assert(Updater<KeepBest>);

} // namespace genetic_algorithm::population_update


#endif // !GA_ALGORITHM_SOGA_UPDATE_HPP