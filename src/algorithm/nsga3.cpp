﻿/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nsga3.hpp"
#include "nd_sort.hpp"
#include "reference_lines.hpp"
#include "../core/ga_info.hpp"
#include "../population/population.hpp"
#include "../metrics/pop_stats.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/math.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <execution>
#include <functional>
#include <iterator>
#include <vector>
#include <span>
#include <utility>
#include <cstddef>

namespace gapp::algorithm
{
    using namespace dtl;
    using math::Point;

    /* Achievement scalarization function. */
    static constexpr double ASF(std::span<const double> ideal_point, std::span<const double> weights, std::span<const double> fitness) noexcept
    {
        GAPP_ASSERT(!ideal_point.empty());
        GAPP_ASSERT(weights.size() == ideal_point.size());
        GAPP_ASSERT(fitness.size() == weights.size());

        double dmax = -math::inf<double>;

        for (size_t i = 0; i < fitness.size(); i++)
        {
            dmax = std::max(dmax, (ideal_point[i] - fitness[i]) / weights[i]);
        }
        
        return dmax;
    }

    /* Create a weight vector for the given axis (used in the ASF). */
    static inline std::vector<double> weightVector(size_t dimensions, size_t axis)
    {
        GAPP_ASSERT(dimensions > axis);

        std::vector weights(dimensions, 1E-6);
        weights[axis] = 1.0;

        return weights;
    }

    /* Find an approximation of the pareto front's nadir point using the minimum of the extreme points. */
    static inline Point findNadirPoint(const std::vector<Point>& extreme_points)
    {
        GAPP_ASSERT(!extreme_points.empty());

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        Point nadir = extreme_points[0];
        for (size_t i = 1; i < extreme_points.size(); i++)
        {
            nadir = detail::elementwise_min(std::move(nadir), extreme_points[i]);
        }

        return nadir;
    }

    /* Normalize a fitness vector using the ideal and nadir points. */
    static inline FitnessVector normalizeFitnessVec(std::span<const double> fvec, std::span<const double> ideal_point, std::span<const double> nadir_point)
    {
        GAPP_ASSERT(fvec.size() == ideal_point.size());
        GAPP_ASSERT(ideal_point.size() == nadir_point.size());

        FitnessVector fnorm(fvec.size());

        for (size_t i = 0; i < fnorm.size(); i++)
        {
            fnorm[i] = (ideal_point[i] - fvec[i]) / std::max(ideal_point[i] - nadir_point[i], 1E-6);
        }

        return fnorm;
    }


    /* NSGA3 IMPLEMENTATION */

    struct NSGA3::Impl
    {
        struct CandidateInfo
        {
            size_t rank;
            size_t ref_idx;
            double ref_dist;
        };

        RefLineGenerator ref_generator_;
        std::vector<Point> ref_lines_;

        std::vector<CandidateInfo> sol_info_;
        std::vector<size_t> niche_counts_;

        Point ideal_point_;
        Point nadir_point_;
        std::vector<Point> extreme_points_;

        /* Generate n reference points in dim dimensions. */
        std::vector<Point> generateReferencePoints(size_t dim, size_t num_points) const;

        /* Update the ideal point approximation using the new points in fmat, assuming maximization. */
        void updateIdealPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Update the extreme points using the new points in fmat, assuming maximization. */
        void updateExtremePoints(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Update the current nadir point based on the extreme points. */
        void updateNadirPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Recalculate the niche counts of the reference lines based on the ref lines associated with the candidates in [first, last). */
        void recalcNicheCounts(ParetoFronts::const_iterator first, ParetoFronts::const_iterator last);

        /* Find the closest reference and its distance for each of the points in the fitness matrix. */
        void associatePopWithRefs(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last,
                                  ParetoFronts::const_iterator pfirst, ParetoFronts::const_iterator plast);

        /* Return true if pop[lhs] is better than pop[rhs]. */
        bool nichedCompare(size_t lhs, size_t rhs) const noexcept;

        /* Return the associated reference direction of a candidate. */
        size_t refIndexOf(const FrontInfo& sol) const noexcept;

        /* Return the associated reference direction's distance for a candidate. */
        double refDistOf(const FrontInfo& sol) const noexcept;

        /* Return the (unique) reference indices that are associated with at least one element in the range [first, last), sorted based on their niche counts. */
        std::vector<size_t> referenceSetOf(ParetoFronts::const_iterator first, ParetoFronts::const_iterator last);

        /* Return the closest solution associated with ref in [first, last). */
        ParetoFronts::iterator findClosestAssociated(ParetoFronts::iterator first, ParetoFronts::iterator last, size_t ref_idx) const noexcept;

        /* Increment the niche count of ref, while keeping refs sorted based on niche counts. */
        void incrementNicheCount(std::vector<size_t>& refs, size_t ref);

        /* Create a new population from the pareto fronts range [first, last). */
        std::vector<size_t> createPopulation(ParetoFronts::const_iterator first, ParetoFronts::const_iterator last);
    };


    NSGA3::NSGA3(RefLineGenerator gen) :
        pimpl_(std::make_unique<Impl>())
    {
        pimpl_->ref_generator_ = std::move(gen);
    }

    NSGA3::NSGA3(const NSGA3& rhs) :
        pimpl_(std::make_unique<Impl>(*rhs.pimpl_))
    {}

    NSGA3& NSGA3::operator=(NSGA3 rhs) noexcept
    {
        this->pimpl_.swap(rhs.pimpl_);
        return *this;
    }

    NSGA3::NSGA3(NSGA3&&) noexcept            = default;
    NSGA3::~NSGA3()                           = default;


    inline std::vector<Point> NSGA3::Impl::generateReferencePoints(size_t dim, size_t num_points) const
    {
        auto ref_lines = ref_generator_(dim, num_points);
        for (auto& ref_line : ref_lines)
        {
            ref_line = math::normalizeVector(std::move(ref_line));
        }

        return ref_lines;
    }

    inline void NSGA3::Impl::updateIdealPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);

        ideal_point_ = detail::elementwise_max(std::move(ideal_point_), detail::maxFitness(first, last));
    }

    void NSGA3::Impl::updateExtremePoints(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);
        GAPP_ASSERT(first->size() == ideal_point_.size());

        const auto popsize = size_t(last - first);

        std::vector<Point> new_extreme_points;
        new_extreme_points.reserve(extreme_points_.size());

        for (size_t dim = 0; dim < ideal_point_.size(); dim++)
        {
            auto weights = weightVector(ideal_point_.size(), dim);
            auto ASFi = [&](const auto& fvec) { return ASF(ideal_point_, weights, fvec); };
            
            std::vector<double> chebysev_distances(popsize + extreme_points_.size());

            std::transform(first, last, chebysev_distances.begin(), ASFi);
            std::transform(extreme_points_.begin(), extreme_points_.end(), chebysev_distances.begin() + popsize, ASFi);

            size_t idx = detail::argmin(chebysev_distances.begin(), chebysev_distances.end());

            (idx >= popsize) ?
                new_extreme_points.push_back(extreme_points_[idx - popsize]) :
                new_extreme_points.push_back(FitnessVector(first[idx]));
        }

        extreme_points_ = std::move(new_extreme_points);
    }

    inline void NSGA3::Impl::updateNadirPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        updateExtremePoints(first, last);
        nadir_point_ = findNadirPoint(extreme_points_);
    }

    inline void NSGA3::Impl::recalcNicheCounts(ParetoFronts::const_iterator first, ParetoFronts::const_iterator last)
    {
        std::fill(niche_counts_.begin(), niche_counts_.end(), 0_sz);
        std::for_each(first, last, [&](const FrontInfo& sol) { niche_counts_[refIndexOf(sol)]++; });
    }

    void NSGA3::Impl::associatePopWithRefs(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last,
                                           ParetoFronts::const_iterator pfirst, ParetoFronts::const_iterator plast)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);
        GAPP_ASSERT(!ref_lines_.empty());

        updateIdealPoint(first, last);
        updateNadirPoint(first, last);

        sol_info_.resize(last - first);

        std::for_each(GAPP_EXEC_UNSEQ, pfirst, plast, [&](const FrontInfo& sol)
        {
            const FitnessVector fnorm = normalizeFitnessVec(first[sol.idx], ideal_point_, nadir_point_);

            auto idistance = [&](const auto& line) { return std::inner_product(fnorm.begin(), fnorm.end(), line.begin(), 0.0); };

            auto closest = detail::max_element(ref_lines_.begin(), ref_lines_.end(), idistance);

            sol_info_[sol.idx].ref_idx  = std::distance(ref_lines_.begin(), closest);
            sol_info_[sol.idx].ref_dist = math::perpendicularDistanceSq(*closest, fnorm);
        });
    }

    inline bool NSGA3::Impl::nichedCompare(size_t lhs, size_t rhs) const noexcept
    {
        /* This version of the compare implementation is from the U-NSGA-III algorithm. */
        if (sol_info_[lhs].ref_idx == sol_info_[rhs].ref_idx)
        {
            if (sol_info_[lhs].rank != sol_info_[rhs].rank)
            {
                return sol_info_[lhs].rank < sol_info_[rhs].rank;
            }

            return sol_info_[lhs].ref_dist < sol_info_[rhs].ref_dist;
        }

        return rng::randomBool();
    }

    inline size_t NSGA3::Impl::refIndexOf(const FrontInfo& sol) const noexcept
    {
        return sol_info_[sol.idx].ref_idx;
    }

    inline double NSGA3::Impl::refDistOf(const FrontInfo& sol) const noexcept
    {
        return sol_info_[sol.idx].ref_dist;
    }

    inline std::vector<size_t> NSGA3::Impl::referenceSetOf(ParetoFronts::const_iterator first, ParetoFronts::const_iterator last)
    {
        std::vector<size_t> refs(last - first);
        std::transform(first, last, refs.begin(), [this](const FrontInfo& sol) { return refIndexOf(sol); });

        detail::erase_duplicates(refs);
        std::sort(refs.begin(), refs.end(), [this](size_t lhs, size_t rhs)
        {
            return niche_counts_[lhs] < niche_counts_[rhs];
        });

        return refs;
    }

    inline ParetoFronts::iterator NSGA3::Impl::findClosestAssociated(ParetoFronts::iterator first, ParetoFronts::iterator last, size_t ref) const noexcept
    {
        auto closest = first;
        auto min_dist = math::inf<double>;

        for (; first != last; ++first)
        {
            if (refIndexOf(*first) == ref && refDistOf(*first) < min_dist)
            {
                closest = first;
                min_dist = refDistOf(*first);
            }
        }

        return closest;
    }

    inline void NSGA3::Impl::incrementNicheCount(std::vector<size_t>& refs, size_t ref)
    {
        niche_counts_[ref]++;

        auto current = std::find(refs.begin(), refs.end(), ref);
        auto first_eq = std::find_if(std::next(current), refs.end(), [&](size_t idx) { return niche_counts_[idx] >= niche_counts_[ref]; });

        std::iter_swap(current, std::prev(first_eq));
    }

    std::vector<size_t> NSGA3::Impl::createPopulation(ParetoFronts::const_iterator first, ParetoFronts::const_iterator last)
    {
        std::vector<size_t> new_pop;
        std::vector<Impl::CandidateInfo> new_info;

        new_pop.reserve(last - first);
        new_info.reserve(last - first);

        for (; first != last; ++first)
        {
            new_pop.push_back(first->idx);
            new_info.push_back(sol_info_[first->idx]);
        }
        sol_info_ = std::move(new_info);

        return new_pop;
    }


    /* NSGA3 INTERFACE FUNCTIONS */

    void NSGA3::initializeImpl(const GaInfo& ga)
    {
        GAPP_ASSERT(ga.population_size() != 0);
        GAPP_ASSERT(ga.num_objectives() > 1, "The number of objectives must be greater than 1 for the NSGA-III algorithm.");

        const auto& fitness_matrix = ga.fitness_matrix();

        pimpl_->ideal_point_ = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());
        pimpl_->extreme_points_ = {};

        pimpl_->ref_lines_ = pimpl_->generateReferencePoints(ga.num_objectives(), ga.population_size());
        pimpl_->niche_counts_.resize(pimpl_->ref_lines_.size());

        auto pfronts = nonDominatedSort(fitness_matrix.begin(), fitness_matrix.end());

        pimpl_->sol_info_.resize(ga.population_size());
        std::for_each(pfronts.begin(), pfronts.end(), [this](const FrontInfo& sol) { pimpl_->sol_info_[sol.idx].rank = sol.rank; });

        pimpl_->associatePopWithRefs(fitness_matrix.begin(), fitness_matrix.end(), pfronts.begin(), pfronts.end());
        pimpl_->recalcNicheCounts(pfronts.begin(), pfronts.end());
    }

    std::vector<size_t> NSGA3::nextPopulationImpl(const GaInfo& ga, FitnessMatrix::const_iterator parents_first,
                                                                    FitnessMatrix::const_iterator /* parents_last */,
                                                                    FitnessMatrix::const_iterator children_last)
    {
        GAPP_ASSERT(ga.num_objectives() > 1);
        GAPP_ASSERT(size_t(children_last - parents_first) >= ga.population_size());
        GAPP_ASSERT(parents_first->size() == ga.num_objectives());

        const size_t popsize = ga.population_size();

        auto pfronts = dtl::nonDominatedSort(parents_first, children_last);
        auto [partial_first, partial_last] = findPartialFront(pfronts.begin(), pfronts.end(), popsize);

        pimpl_->sol_info_.resize(children_last - parents_first);
        std::for_each(pfronts.begin(), pfronts.end(), [this](const FrontInfo& sol) { pimpl_->sol_info_[sol.idx].rank = sol.rank; });

        /* The ref lines of the candidates after partial_last are irrelevant, as they can never be part of the next population. */
        pimpl_->associatePopWithRefs(parents_first, children_last, pfronts.begin(), partial_last);
        /* The niche counts should be calculated excluding the partial front for now. */
        pimpl_->recalcNicheCounts(pfronts.begin(), partial_first);

        /* Find the reference lines associated with the partial front. */
        std::vector<size_t> refs = pimpl_->referenceSetOf(partial_first, partial_last);

        /* Move the best elements to the front of the partial front. */
        while (partial_first != pfronts.begin() + popsize)
        {
            const size_t min_count = pimpl_->niche_counts_[refs[0]];
            const auto min_last = std::find_if(refs.begin(), refs.end(), [&](size_t idx) { return pimpl_->niche_counts_[idx] != min_count; });

            size_t ref = *rng::randomElement(refs.begin(), min_last);
            pimpl_->incrementNicheCount(refs, ref);

            const size_t assoc_count = std::count_if(partial_first, partial_last, [&](const FrontInfo& sol) { return pimpl_->refIndexOf(sol) == ref; });
            const auto closest = pimpl_->findClosestAssociated(partial_first, partial_last, ref);

            /* Move the selected candidate to the front of the partial front so it can't be selected again. */
            std::iter_swap(closest, partial_first++);

            /* If the selected candidate was the only one in the partial front associated with this reference direction,
               the reference direction needs to also be removed. */
            if (assoc_count == 1) detail::erase_first_stable(refs, ref);
        }

        return pimpl_->createPopulation(pfronts.begin(), pfronts.begin() + popsize);
    }

    size_t NSGA3::selectImpl(const GaInfo&, const FitnessMatrix& fmat) const
    {
        GAPP_ASSERT(!fmat.empty());

        const size_t idx1 = rng::randomIdx(fmat);
        const size_t idx2 = rng::randomIdx(fmat);

        return pimpl_->nichedCompare(idx1, idx2) ? idx1 : idx2;
    }

    std::vector<size_t> NSGA3::optimalSolutionsImpl(const GaInfo&) const
    {
        return detail::find_indices(pimpl_->sol_info_, [](const Impl::CandidateInfo& sol) { return sol.rank == 0; });
    }

} // namespace gapp::algorithm