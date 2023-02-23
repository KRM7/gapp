/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_NSGA3_HPP
#define GA_ALGORITHM_NSGA3_HPP

#include "algorithm_base.hpp"
#include "reference_lines.hpp"
#include "../utility/math.hpp"
#include <vector>
#include <functional>
#include <memory>
#include <cstddef>

namespace genetic_algorithm::algorithm
{
    /**
    * NSGA-III algorithm, used for many-objective optimization. (Doesn't work for single-objective problems.) \n
    * The aim of the algorithm is to find a set of solutions which are well-spread out
    * along the entire pareto-front (in the objective-space). \n
    * 
    * The algorithm uses non-dominated sorting to sort the solutions into pareto fronts,
    * and then selects the candidates of the best fronts for the population of the next generation. \n
    * Candidates that belong to the same front are ranked using a set of reference
    * directions in the objective space. Candidates closest to reference directions which have less candidates associated with them,
    * and candidates closer to reference directions are considered better. \n
    * The reference directions are generated at the start of the run, and don't change throughout it. \n
    * 
    * The algorithm uses a selection operator that selects solutions based on these same criteria
    * (their pareto ranks and their distances from the reference directions). \n
    * 
    * The selection and population update methods of this algorithm can't be changed. \n
    * Has no parameters. \n
    * 
    * See:
    *  Deb, K., and Jain, H. "An evolutionary many-objective optimization algorithm using reference-point-based nondominated sorting approach, part I:
    *  solving problems with box constraints." IEEE transactions on evolutionary computation 18, no. 4 (2013): 577-601.
    */
    class NSGA3 final : public Algorithm
    {
    public:
        /** The type of the reference line generator function. */
        using RefLineGenerator = std::function<std::vector<math::Point>(size_t, size_t)>;

        /**
        * Create an NSGA-III algorithm.
        * 
        * @param gen The method to use for generating the reference lines for the algorithm.
        */
        NSGA3(RefLineGenerator gen = reflines::quasirandomSimplexPointsMirror);

        NSGA3(const NSGA3&);
        NSGA3(NSGA3&&) noexcept;
        NSGA3& operator=(NSGA3) noexcept;
        ~NSGA3() override;

    private:

        void initializeImpl(const GaInfo& ga) override;
        size_t selectImpl(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<size_t> nextPopulationImpl(const GaInfo& ga,
                                               FitnessMatrix::const_iterator first,
                                               FitnessMatrix::const_iterator children_first,
                                               FitnessMatrix::const_iterator last) override;

        std::vector<size_t> optimalSolutionsImpl(const GaInfo& ga) const override;

        struct Impl;
        std::unique_ptr<Impl> pimpl_;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_NSGA2_HPP