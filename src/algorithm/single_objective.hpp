/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_SINGLE_OBJECTIVE_HPP
#define GA_ALGORITHM_SINGLE_OBJECTIVE_HPP

#include "algorithm_base.decl.hpp"
#include "selection_base.hpp"
#include "replacement_base.hpp"
#include "soga_selection.hpp"
#include "soga_replacement.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <concepts>
#include <functional>
#include <memory>
#include <cstddef>

namespace gapp::algorithm
{
    /**
    * A generic algorithm for single-objective optimization.
    * 
    * The algorithm combines a selection method and a population replacement
    * method. The selection method is used to select candidates from the populations
    * for crossover, while the population replacement method is used to create the population
    * for the next generation of the algorithm from the combined parent and child populations.
    * 
    * Move-only.
    */
    class SingleObjective final : public Algorithm
    {
    public:
        using DefaultSelection   = selection::Tournament; /**< The selection method used when not specified explicitly. */
        using DefaultReplacement = replacement::KeepBest; /**< The population update method used when not specified explicitly. */

        /**
        * The general callable type that can be used as a selection method,
        * when not using a selection method derived from selection::Selection.
        * @see selection_method()
        */
        using SelectionCallable = std::function<size_t(const GaInfo&, const FitnessMatrix&)>;

        /**
        * The general callable type that can be used as a population replacement policy,
        * when not using a replacement policy derived from replacement::Replacement.
        * @see replacement_method()
        */
        using ReplacementCallable = std::function<std::vector<size_t>(const GaInfo&, const FitnessMatrix&)>;


        /**
        * Create a single-objective algorithm.
        *
        * @param selection The selection method to use.
        * @param replacement The replacement policy to use.
        */
        template<typename S = DefaultSelection, typename R = DefaultReplacement>
        requires std::derived_from<S, selection::Selection> && std::derived_from<R, replacement::Replacement>
        explicit SingleObjective(S selection = S{}, R replacement = R{});

        /**
        * Create a single-objective algorithm.
        *
        * @param selection The selection method to use. Can't be a nullptr.
        * @param replacement The replacement policy to use. Can't be a nullptr.
        */
        explicit SingleObjective(std::unique_ptr<selection::Selection> selection, std::unique_ptr<replacement::Replacement> replacement = std::make_unique<DefaultReplacement>());

        /**
        * Create a single-objective algorithm using the default replacement method.
        *
        * @param selection The selection method to use. Can't be a nullptr.
        */
        explicit SingleObjective(SelectionCallable selection);

        /**
        * Create a single-objective algorithm.
        *
        * @param selection The selection method to use. Can't be a nullptr.
        * @param replacement The replacement policy to use. Can't be a nullptr.
        */
        SingleObjective(SelectionCallable selection, ReplacementCallable replacement);

        
        /**
        * Set the selection method used by the algorithm.
        *
        * @param selection The selection method to use.
        */
        template<typename S>
        requires std::derived_from<S, selection::Selection>
        void selection_method(S selection);

        /**
        * Set the selection method used by the algorithm.
        *
        * @param selection The selection method to use. Can't be a nullptr.
        */
        void selection_method(std::unique_ptr<selection::Selection> selection);

        /**
        * Set the selection method used by the algorithm.
        * The function used should be thread-safe if parallel execution is enabled
        * (which is true by default).
        *
        * @param f The selection method to use. Can't be a nullptr.
        */
        void selection_method(SelectionCallable f);

        /** @returns The selection operator used by the algorithm. */
        [[nodiscard]]
        selection::Selection& selection_method() & noexcept { GAPP_ASSERT(selection_); return *selection_; }

        /** @returns The selection operator used by the algorithm. */
        [[nodiscard]]
        const selection::Selection& selection_method() const& noexcept { GAPP_ASSERT(selection_); return *selection_; }


        /**
        * Set the population replacement policy used by the algorithm.
        *
        * @param replacement The replacement policy to use.
        */
        template<typename R>
        requires std::derived_from<R, replacement::Replacement>
        void replacement_method(R replacement);

        /**
        * Set the population replacement policy used by the algorithm.
        *
        * @param replacement The replacement policy to use. Can't be a nullptr.
        */
        void replacement_method(std::unique_ptr<replacement::Replacement> replacement);

        /**
        * Set the population replacement policy used by the algorithm.
        *
        * @param f The replacement policy to use. Can't be a nullptr.
        */
        void replacement_method(ReplacementCallable f);

        /** @returns The population replacement policy used by the algorithm. */
        [[nodiscard]]
        replacement::Replacement& replacement_method() & noexcept { GAPP_ASSERT(replacement_); return *replacement_; }

        /** @returns The population replacement policy used by the algorithm. */
        [[nodiscard]]
        const replacement::Replacement& replacement_method() const& noexcept { GAPP_ASSERT(replacement_); return *replacement_; }

    private:

        void initializeImpl(const GaInfo& ga) override;
        void prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<size_t> nextPopulationImpl(const GaInfo& ga, const FitnessMatrix& fmat) override;

        std::unique_ptr<selection::Selection> selection_;
        std::unique_ptr<replacement::Replacement> replacement_;
    };

} // namespace gapp::algorithm


/* IMPLEMENTATION */

#include <utility>

namespace gapp::algorithm
{
    template<typename S, typename R>
    requires std::derived_from<S, selection::Selection> && std::derived_from<R, replacement::Replacement>
    inline SingleObjective::SingleObjective(S selection, R replacement) :
        selection_(std::make_unique<S>(std::move(selection))), replacement_(std::make_unique<R>(std::move(replacement)))
    {}

    template<typename S>
    requires std::derived_from<S, selection::Selection>
    inline void SingleObjective::selection_method(S selection)
    {
        selection_ = std::make_unique<S>(std::move(selection));
    }

    template<typename R>
    requires std::derived_from<R, replacement::Replacement>
    inline void SingleObjective::replacement_method(R replacement)
    {
        replacement_ = std::make_unique<R>(std::move(replacement));
    }

    inline void SingleObjective::prepareSelectionsImpl(const GaInfo& ga, const FitnessMatrix& fmat)
    {
        GAPP_ASSERT(selection_);

        selection_->prepareSelectionsImpl(ga, fmat);
    }

    inline size_t SingleObjective::selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const
    {
        GAPP_ASSERT(selection_);

        return selection_->selectImpl(ga, fmat);
    }

} // namespace gapp::algorithm

#endif // !GA_ALGORITHM_SINGLE_OBJECTIVE_HPP