/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ALGORITHM_NSGA3_HPP
#define GA_ALGORITHM_NSGA3_HPP

#include "algorithm_base.hpp"

namespace genetic_algorithm::algorithm
{
    /**
    * NSGA-III algorithm for many-objective optimization.
    */
    class NSGA3 final : public Algorithm
    {
    public:
        void initialize(const GaInfo& ga) override;
        void prepareSelections(const GaInfo&, const FitnessMatrix&) override {}
        size_t select(const GaInfo& ga, const FitnessMatrix& fmat) override;
        std::vector<size_t> nextPopulation(const GaInfo& ga,
                                           FitnessMatrix::const_iterator first,
                                           FitnessMatrix::const_iterator children_first,
                                           FitnessMatrix::const_iterator last) override;

        using Point = std::vector<double>;

        struct RefPoint
        {
            const Point point;
            size_t niche_count = 0;

            RefPoint(const Point& p) : point(p) {}
        };

    private:

        struct CandidateInfo
        {
            size_t rank = 0;
            size_t ref_idx = 0;
            double ref_dist = 0.0;
        };

        using ASF = std::function<double(const std::vector<double>&)>;

        std::vector<CandidateInfo> sol_info_;
        std::vector<RefPoint> ref_points_;

        Point ideal_point_;
        Point nadir_point_;
        std::vector<Point> extreme_points_;

        /* Generate n reference points on the unit simplex in dim dimensions from a uniform distribution. */
        static std::vector<Point> generateRefPoints(size_t n, size_t dim);

        /* Find the index and distance of the closest reference line to the point p. */
        static std::pair<size_t, double> findClosestRef(const std::vector<RefPoint>& refs, const Point& p);

        /* Create an achievement scalarization function. */
        static ASF getASF(std::vector<double> ideal_point, std::vector<double> weights) noexcept;

        /* Create a weight vector for the given axis (used in the ASF). */
        static std::vector<double> weightVector(size_t dimensions, size_t axis);

        /* Update the approximate ideal point using new points in fmat, assuming maximization. */
        static void updateIdealPoint(Point& ideal_point, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Update the extreme points with the new points in fmat. */
        static void updateExtremePoints(std::vector<Point>& extreme_points, const Point& ideal_point, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Find an approximation of the nadir point of the pareto front using the minimum of the extreme points. */
        static Point findNadirPoint(const std::vector<Point>& extreme_points);

        /* TODO ... map fitness values onto unit simplex */
        static FitnessVector normalize(const FitnessVector& fvec, const Point& ideal_point, const Point& nadir_point);

        /* Find the closest reference point to each candidate after normalization, and their distances. */
        void associatePopWithRefs(std::vector<CandidateInfo>& props, FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, const std::vector<RefPoint>& refs);

        /* Returns true if Pop[lhs] is better than Pop[rhs]. */
        bool nichedCompare(size_t lhs, size_t rhs) const noexcept;

        /* Return the niche counts of the ref points and assign niche counts to the candidates. */
        static void updateNicheCounts(std::vector<RefPoint>& refs, const std::vector<CandidateInfo>& props) noexcept;

        /* Returns the niche counts of the given candidate. */
        size_t& nicheCountOf(const CandidateInfo& info) noexcept;
        size_t& nicheCountOf(size_t sol_idx) noexcept;
        const size_t& nicheCountOf(const CandidateInfo& info) const noexcept;
        const size_t& nicheCountOf(size_t sol_idx) const noexcept;
    };

} // namespace genetic_algorithm::algorithm

#endif // !GA_ALGORITHM_NSGA2_HPP