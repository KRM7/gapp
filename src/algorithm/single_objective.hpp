/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_HPP

#include "algorithm_base.decl.hpp"
#include "soga_selection.hpp"
#include "soga_update.hpp"
#include <functional>
#include <memory>
#include <functional>
#include <type_traits>
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    /**
    * A generic algorithm for single-objective optimization. \n
    * The algorithm combines a selection method and a population update
    * method. The selection method is used to select candidates from the populations
    * for crossover, while the population update method is used to create the population
    * for the next generation of the algorithm from the combined parent and child populations.
    * 
    * Move-only.
    */
    class SingleObjective final : public Algorithm
    {
    public:
        using DefaultSelection = selection::Tournament; /**< The selection method used when not specified explicitly. */
        using DefaultUpdater   = update::KeepBest;      /**< The population method used when not specified explicitly. */

        /** The signature of the selection function. */
        using SelectionFunction = std::function<size_t(const GaInfo&, const FitnessMatrix&)>;

        /** The signature of the population update function. */
        using UpdateFunction = std::function<std::vector<size_t>(const GaInfo&, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator)>;

        /**
        * Create a single objective algorithm using the default selection
        * and population update methods.
        */
        SingleObjective();

        /**
        * Create a single objective algorithm using the default population update method.
        *
        * @param selection The selection method to use in the algorithm.
        */
        template<typename S>
        requires selection::SelectionType<S> && std::is_final_v<S>
        explicit SingleObjective(S selection);

        /**
        * Create a single objective algorithm.
        *
        * @param selection The selection method to use in the algorithm.
        * @param updater The method used to update the population between generations of the algorithm.
        */
        template<typename S, typename U>
        requires selection::SelectionType<S> && std::is_final_v<S> && update::UpdaterType<U> && std::is_final_v<U>
        SingleObjective(S selection, U updater);

        /**
        * Create a single objective algorithm using the default population update method.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        */
        template<selection::SelectionType S>
        explicit SingleObjective(std::unique_ptr<S>&& selection);

        /**
        * Create a single objective algorithm.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        * @param updater The method used to update the population between generations of the algorithm. Can't be a nullptr.
        */
        template<selection::SelectionType S, update::UpdaterType U>
        SingleObjective(std::unique_ptr<S>&& selection, std::unique_ptr<U>&& updater);

        /**
        * Create a single objective algorithm using the default population update method.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        */
        explicit SingleObjective(SelectionFunction selection);

        /**
        * Create a single objective algorithm.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        * @param updater The method used to update the population between generations of the algorithm. Can't be a nullptr.
        */
        SingleObjective(SelectionFunction selection, UpdateFunction updater);

        
        /**
        * Set the selection method used by the algorithm.
        * 
        * @see Selection
        *
        * @param selection The selection method used by the algorithm.
        */
        template<typename S>
        requires selection::SelectionType<S> && std::is_final_v<S>
        void selection_method(S selection);

        /**
        * Set the selection method used by the algorithm.
        * 
        * @see Selection
        *
        * @param selection The selection method used by the algorithm. Can't be a nullptr.
        */
        template<selection::SelectionType S>
        void selection_method(std::unique_ptr<S>&& selection);

        /**
        * Set the selection method used by the algorithm. \n
        * The function used should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Selection
        * @see SelectionFunction
        *
        * @param f The selection method used by the algorithm. Can't be a nullptr.
        */
        void selection_method(SelectionFunction f);

        /** @returns The selection operator used by the algorithm. */
        [[nodiscard]]
        selection::Selection& selection_method() & noexcept { return *selection_; }

        /** @returns The selection operator used by the algorithm. */
        [[nodiscard]]
        const selection::Selection& selection_method() const& noexcept { return *selection_; }


        /**
        * Set the population update method used by the algorithm.
        * 
        * @see Updater
        *
        * @param updater The population update method used by the algorithm.
        */
        template<typename U>
        requires update::UpdaterType<U> && std::is_final_v<U>
        void update_method(U updater);

        /**
        * Set the population update method used by the algorithm.
        * 
        * @see Updater
        *
        * @param updater The population update method used by the algorithm. Can't be a nullptr.
        */
        template<update::UpdaterType U>
        void update_method(std::unique_ptr<U>&& updater);

        /**
        * Set the population update method used by the algorithm.
        *
        * @see Updater
        * @see UpdateFunction
        *
        * @param f The population update method used by the algorithm. Can't be a nullptr.
        */
        void update_method(UpdateFunction f);

        /** @returns The population update operator used by the algorithm. */
        [[nodiscard]]
        update::Updater& update_method() & noexcept { return *updater_; }

        /** @returns The population update operator used by the algorithm. */
        [[nodiscard]]
        const update::Updater& update_method() const& noexcept { return *updater_; }

    private:

        void initializeImpl(const GaInfo& ga) override;
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator children_first, FitnessMatrix::const_iterator last) override;

        std::unique_ptr<selection::Selection> selection_;
        std::unique_ptr<update::Updater> updater_;
    };

} // namespace genetic_algorithm::algorithm


/* IMPLEMENTATION */

#include "../utility/utility.hpp"
#include <stdexcept>
#include <utility>

namespace genetic_algorithm::algorithm
{
    inline SingleObjective::SingleObjective() :
        Algorithm(), selection_(std::make_unique<DefaultSelection>()), updater_(std::make_unique<DefaultUpdater>())
    {}

    template<typename S>
    requires selection::SelectionType<S> && std::is_final_v<S>
    inline SingleObjective::SingleObjective(S selection) :
        SingleObjective(std::move(selection), DefaultUpdater{})
    {}

    template<typename S, typename U>
    requires selection::SelectionType<S> && std::is_final_v<S> && update::UpdaterType<U> && std::is_final_v<U>
    inline SingleObjective::SingleObjective(S selection, U updater) :
        Algorithm(), selection_(std::make_unique<S>(std::move(selection))), updater_(std::make_unique<U>(std::move(updater)))
    {}

    template<selection::SelectionType S>
    inline SingleObjective::SingleObjective(std::unique_ptr<S>&& selection) :
        SingleObjective(std::move(selection), std::make_unique<DefaultUpdater>())
    {}

    template<selection::SelectionType S, update::UpdaterType U>
    inline SingleObjective::SingleObjective(std::unique_ptr<S>&& selection, std::unique_ptr<U>&& updater) :
        Algorithm(), selection_(std::move(selection)), updater_(std::move(updater))
    {
        if (!selection_) GA_THROW(std::invalid_argument, "The selection method can't be a nullptr.");
        if (!updater_) GA_THROW(std::invalid_argument, "The population update method can't be a nullptr.");
    }

    template<typename S>
    requires selection::SelectionType<S> && std::is_final_v<S>
    inline void SingleObjective::selection_method(S selection)
    {
        selection_ = std::make_unique<S>(std::move(selection));
    }

    template<selection::SelectionType S>
    inline void SingleObjective::selection_method(std::unique_ptr<S>&& selection)
    {
        if (!selection_) GA_THROW(std::invalid_argument, "The selection method can't be a nullptr.");

        selection_ = std::move(selection);
    }

    template<typename U>
    requires update::UpdaterType<U> && std::is_final_v<U>
    inline void SingleObjective::update_method(U updater)
    {
        updater_ = std::make_unique<U>(std::move(updater));
    }

    template<update::UpdaterType U>
    inline void SingleObjective::update_method(std::unique_ptr<U>&& updater)
    {
        if (!updater_) GA_THROW(std::invalid_argument, "The population update method can't be a nullptr.");

        updater_ = std::move(updater);
    }

    inline void SingleObjective::prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        selection_->prepareSelectionsImpl(ga, fmat);
    }

    inline size_t SingleObjective::selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const
    {
        return selection_->selectImpl(ga, fmat);
    }

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_HPP