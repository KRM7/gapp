/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_MULTI_OBJECTIVE_HPP
#define GA_SELECTION_MULTI_OBJECTIVE_HPP

#include "selection_base.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::selection::multi_objective
{
    /**
    * Non-dominated sorting genetic algorithm (NSGA-II) for multi-objective optimization.
    */
    class NSGA2 final : public Selection
    {
    public:
        void init(const GaInfo& ga) override;
        size_t select(const GaInfo& ga, const FitnessMatrix& pop) override;
        std::vector<size_t> nextPopulation(const GaInfo& ga, FitnessMatrix& combined_pop) override;

    private:
        std::vector<size_t> ranks_;
        std::vector<double> dists_;

        /* Returns true if Pop[lidx] is better than Pop[ridx]. */
        bool crowdedCompare(size_t lidx, size_t ridx) const noexcept;
    };

    /**
    * NSGA-III algorithm for many-objective optimization.
    */
    class NSGA3 final : public Selection
    {
    public:
        void init(const GaInfo& ga) override;
        size_t select(const GaInfo& ga, const FitnessMatrix& pop) override;
        std::vector<size_t> nextPopulation(const GaInfo& ga, FitnessMatrix& combined_pop) override;

        using Point = std::vector<double>;

        struct RefPoint
        {
            const Point point;
            size_t niche_count;

            RefPoint(const Point& p) : point(p), niche_count(0) {}
        };

    private:

        struct CandidateInfo
        {
            size_t rank = 0;
            size_t ref_idx = 0;
            double ref_dist = 0.0;
        };

        std::vector<CandidateInfo> sol_info_;
        std::vector<RefPoint> ref_points_;

        Point ideal_point_;
        Point nadir_point_;
        std::vector<Point> extreme_points_;

        /* Update the approximate ideal point using new points in fmat, assuming maximization. */
        static void updateIdealPoint(Point& ideal_point, const FitnessMatrix& fmat);

        /* Create a weight vector for the given axis (used in the ASF). */
        static std::vector<double> weightVector(size_t dimensions, size_t axis);

        /* Update the extreme points with the new points in fmat. */
        static void updateExtremePoints(std::vector<Point>& extreme_points, const FitnessMatrix& pop, const Point& ideal_point);

        /* Find an approximation of the nadir point of the pareto front using the minimum of the extreme points. */
        static Point findNadirPoint(const std::vector<Point>& extreme_points);

        /* Returns true if Pop[lhs] is better than Pop[rhs]. */
        bool nichedCompare(size_t lhs, size_t rhs) const noexcept;

        static FitnessVector normalize(const FitnessVector& fvec, const Point& ideal_point, const Point& nadir_point);

        /* Find the closest reference point to each candidate after normalization, and their distances. */
        void associatePopWithRefs(std::vector<CandidateInfo>& props, const FitnessMatrix& fmat, const std::vector<RefPoint>& refs);
        
        /* Returns the niche counts of the given candidate. */
        size_t& nicheCountOf(CandidateInfo& info);
        size_t& nicheCountOf(size_t sol_idx);
        const size_t& nicheCountOf(const CandidateInfo& info) const;
        const size_t& nicheCountOf(size_t sol_idx) const;

        /* Return the niche counts of the ref points and assign niche counts to the candidates. */
        static void updateNicheCounts(std::vector<RefPoint>& refs, const std::vector<CandidateInfo>& props);
    };

} // namespace genetic_algorithm::selection::multi_objective

#endif // !GA_SELECTION_MULTI_OBJECTIVE_HPP