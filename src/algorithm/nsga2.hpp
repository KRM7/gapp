/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_NSGA2_HPP
#define GA_ALGORITHM_NSGA2_HPP

#include "algorithm_base.hpp"

namespace genetic_algorithm::algorithm
{
    /**
    * Non-dominated sorting genetic algorithm (NSGA-II) for multi-objective optimization.
    */
    class NSGA2 final : public Algorithm
    {
    public:
        void initialize(const GaInfo& ga) override;
        void prepareSelections(const GaInfo&, const FitnessMatrix&) override {}
        size_t select(const GaInfo& ga, const FitnessMatrix& pop) override;
        std::vector<size_t> nextPopulation(const GaInfo& ga,
                                           FitnessMatrix::const_iterator first,
                                           FitnessMatrix::const_iterator children_first,
                                           FitnessMatrix::const_iterator last) override;

    private:
        std::vector<size_t> ranks_;
        std::vector<double> dists_;

        /* Returns true if Pop[lidx] is better than Pop[ridx]. */
        bool crowdedCompare(size_t lidx, size_t ridx) const noexcept;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_NSGA2_HPP