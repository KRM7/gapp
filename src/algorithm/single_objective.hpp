/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_HPP

#include "algorithm_base.decl.hpp"
#include "selection_base.hpp"
#include "updater_base.hpp"
#include "soga_selection.hpp"
#include "soga_update.hpp"
#include <vector>
#include <functional>
#include <memory>
#include <cstddef>
#include <cassert>

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
        using DefaultUpdater   = update::KeepBest;      /**< The population update method used when not specified explicitly. */

        /** The general callable type that can be used as a selection method. */
        using SelectionCallable = std::function<size_t(const GaInfo&, const FitnessMatrix&)>;

        /** The general callable type that can be used as a population update method. */
        using UpdateCallable = std::function<std::vector<size_t>(const GaInfo&, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator, FitnessMatrix::const_iterator)>;

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
        explicit SingleObjective(std::unique_ptr<selection::Selection> selection);

        /**
        * Create a single objective algorithm.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        * @param updater The method used to update the population between generations of the algorithm. Can't be a nullptr.
        */
        SingleObjective(std::unique_ptr<selection::Selection> selection, std::unique_ptr<update::Updater> updater);

        /**
        * Create a single objective algorithm using the default population update method.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        */
        explicit SingleObjective(SelectionCallable selection);

        /**
        * Create a single objective algorithm.
        *
        * @param selection The selection method to use in the algorithm. Can't be a nullptr.
        * @param updater The method used to update the population between generations of the algorithm. Can't be a nullptr.
        */
        SingleObjective(SelectionCallable selection, UpdateCallable updater);

        
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
        void selection_method(std::unique_ptr<selection::Selection> selection);

        /**
        * Set the selection method used by the algorithm. \n
        * The function used should be thread-safe if parallel execution is enabled (enabled by default).
        * 
        * @see Selection
        * @see SelectionCallable
        *
        * @param f The selection method used by the algorithm. Can't be a nullptr.
        */
        void selection_method(SelectionCallable f);

        /** @returns The selection operator used by the algorithm. */
        [[nodiscard]]
        selection::Selection& selection_method() & noexcept { assert(selection_); return *selection_; }

        /** @returns The selection operator used by the algorithm. */
        [[nodiscard]]
        const selection::Selection& selection_method() const& noexcept { assert(selection_); return *selection_; }


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
        void update_method(std::unique_ptr<update::Updater> updater);

        /**
        * Set the population update method used by the algorithm.
        *
        * @see Updater
        * @see UpdateCallable
        *
        * @param f The population update method used by the algorithm. Can't be a nullptr.
        */
        void update_method(UpdateCallable f);

        /** @returns The population update operator used by the algorithm. */
        [[nodiscard]]
        update::Updater& update_method() & noexcept { assert(updater_); return *updater_; }

        /** @returns The population update operator used by the algorithm. */
        [[nodiscard]]
        const update::Updater& update_method() const& noexcept { assert(updater_); return *updater_; }

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

#include <utility>

namespace genetic_algorithm::algorithm
{
    inline SingleObjective::SingleObjective() :
        selection_(std::make_unique<DefaultSelection>()), updater_(std::make_unique<DefaultUpdater>())
    {}

    template<typename S>
    requires selection::SelectionType<S> && std::is_final_v<S>
    inline SingleObjective::SingleObjective(S selection) :
        SingleObjective(std::move(selection), DefaultUpdater{})
    {}

    template<typename S, typename U>
    requires selection::SelectionType<S> && std::is_final_v<S> && update::UpdaterType<U> && std::is_final_v<U>
    inline SingleObjective::SingleObjective(S selection, U updater) :
        selection_(std::make_unique<S>(std::move(selection))), updater_(std::make_unique<U>(std::move(updater)))
    {}

    template<typename S>
    requires selection::SelectionType<S> && std::is_final_v<S>
    inline void SingleObjective::selection_method(S selection)
    {
        selection_ = std::make_unique<S>(std::move(selection));
    }

    template<typename U>
    requires update::UpdaterType<U> && std::is_final_v<U>
    inline void SingleObjective::update_method(U updater)
    {
        updater_ = std::make_unique<U>(std::move(updater));
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