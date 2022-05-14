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
        void prepare(const GaInfo& ga, const FitnessMatrix& pop) override;
        size_t select(const GaInfo& ga, const FitnessMatrix& pop) override;
        std::vector<size_t> nextPopulation(const GaInfo& ga, FitnessMatrix& combined_pop) override;

    private:
        std::vector<size_t> ranks_;
        std::vector<double> dists_;

        /* Returns true if Pop[lidx] is better than Pop[ridx]. */
        constexpr bool crowdedCompare(size_t lidx, size_t ridx) const noexcept;
    };

    /**
    * NSGA-III algorithm for many-objective optimization.
    */
    class NSGA3 final : public Selection
    {
    public:
        void init(const GaInfo& ga) override;
        void prepare(const GaInfo& ga, const FitnessMatrix& pop) override;
        size_t select(const GaInfo& ga, const FitnessMatrix& pop) override;
        std::vector<size_t> nextPopulation(const GaInfo& ga, FitnessMatrix& combined_pop) override;

    private:
        using Point = std::vector<double>;

        struct CandidateInfo
        {
            size_t rank = 0;
            size_t ref_idx = 0;
            double ref_dist = 0.0;
            size_t niche_count = 0;
        };

        std::vector<CandidateInfo> sol_props_;
        std::vector<Point> ref_points_;
        std::vector<size_t> ref_niche_counts_;

        Point ideal_point_;
        Point nadir_point_;
        std::vector<Point> extreme_points_;

        static void updateIdealPoint(Point& ideal_point, const FitnessMatrix& pop) noexcept;
        static std::vector<Point> initExtremePoints(const FitnessMatrix& pop, const Point& ideal_point);
        static void updateExtremePoints(std::vector<Point>& extreme_points, const FitnessMatrix& pop, const Point& ideal_point);
        static Point nadirPoint(const std::vector<Point>& extreme_points);

        /* Find the closest reference point to each candidate after normalization, and their distances. */
        void associatePopWithRefs(const FitnessMatrix& pop);

        /* Return the niche counts of the ref points and assign niche counts to the candidates. */
        static std::vector<size_t> calcNicheCounts(const GaInfo& ga, std::vector<CandidateInfo>& props);

        /* Returns true if Pop[lidx] is better than Pop[ridx]. */
        constexpr bool nichedCompare(size_t lidx, size_t ridx) const noexcept;
    };

} // namespace genetic_algorithm::selection::multi_objective

#endif // !GA_SELECTION_MULTI_OBJECTIVE_HPP