/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_NSGA3_HPP
#define GA_ALGORITHM_NSGA3_HPP

#include "algorithm_base.hpp"

namespace genetic_algorithm::algorithm
{
    namespace dtl { struct FrontInfo; }

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
        void initialize(const GaInfo& ga) override;

        void prepareSelections(const GaInfo&, const FitnessMatrix&) override {}

        size_t select(const GaInfo& ga, const FitnessMatrix& fmat) const override;

        std::vector<size_t> nextPopulation(const GaInfo& ga,
                                           FitnessMatrix::const_iterator first,
                                           FitnessMatrix::const_iterator children_first,
                                           FitnessMatrix::const_iterator last) override;

        std::optional<std::vector<size_t>> optimalSolutions(const GaInfo& ga) const override;

        /* A reference point/direction with its associated niche-count. */
        struct RefPoint
        {
            const std::vector<double> point;
            size_t niche_count = 0;

            RefPoint(std::vector<double> p) : point(std::move(p)) {}
        };

    private:
        /* Stats associated with each of the solutions. */
        struct CandidateInfo
        {
            size_t rank = 0;
            size_t ref_idx = 0;
            double ref_dist = 0.0;
        };

        using ASF   = std::function<double(const std::vector<double>&)>;  /* Achievement scalarization function type. */
        using Point = std::vector<double>;

        std::vector<CandidateInfo> sol_info_;
        std::vector<RefPoint> ref_points_;

        Point ideal_point_;
        Point nadir_point_;
        std::vector<Point> extreme_points_;

        /* Generate n reference points on the unit simplex in dim dimensions from a uniform distribution. */
        static std::vector<RefPoint> generateRefPoints(size_t n, size_t dim);

        /* Find the index and distance of the closest reference line to the point p. */
        static std::pair<size_t, double> findClosestRef(const std::vector<RefPoint>& refs, const Point& p);

        /* Create an achievement scalarization function. */
        static ASF getASF(std::vector<double> ideal_point, std::vector<double> weights) noexcept;

        /* Create a weight vector for the given axis (used in the ASF). */
        static std::vector<double> weightVector(size_t dimensions, size_t axis);

        /* Update the ideal point approximation using the new points in fmat, assuming maximization. */
        void updateIdealPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Update the extreme points using the new points in fmat, assuming maximization. */
        void updateExtremePoints(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Find an approximation of the pareto front's nadir point using the minimum of the extreme points. */
        Point findNadirPoint(const std::vector<Point>& extreme_points);

        /* Normalize a fitness vector using the ideal and nadir points. */
        FitnessVector normalize(const FitnessVector& fvec);

        /* Find the closest reference and its distance for each of the points in the fitness matrix. */
        void associatePopWithRefs(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Return true if pop[lhs] is better than pop[rhs]. */
        bool nichedCompare(size_t lhs, size_t rhs) const noexcept;

        /* Return the associated reference point of a candidate. */
        RefPoint& refPointOf(const dtl::FrontInfo& sol) noexcept;
        const RefPoint& refPointOf(const dtl::FrontInfo& sol) const noexcept;

        /* Return the associated reference point's distance for a candidate. */
        double refDistOf(const dtl::FrontInfo& sol) const noexcept;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_NSGA2_HPP