/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "nsga3.hpp"
#include "nd_sort.hpp"
#include "reference_lines.hpp"
#include "../core/ga_info.hpp"
#include "../core/population.hpp"
#include "../metrics/pop_stats.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include "../utility/small_vector.hpp"
#include "../utility/thread_pool.hpp"
#include "../utility/math.hpp"
#include "../utility/rng.hpp"
#include "../utility/utility.hpp"
#include <algorithm>
#include <functional>
#include <iterator>
#include <vector>
#include <span>
#include <utility>
#include <cstddef>

namespace gapp::algorithm
{
    using namespace dtl;

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
    static inline small_vector<double> weightVector(size_t dimensions, size_t axis)
    {
        GAPP_ASSERT(dimensions > axis);

        small_vector weights(dimensions, 1E-6);
        weights[axis] = 1.0;

        return weights;
    }

    /* Normalize a fitness vector using the ideal and nadir points. */
    static inline FitnessVector normalizeFitnessVec(std::span<const double> fvec, std::span<const double> ideal_point, std::span<const double> nadir_point)
    {
        GAPP_ASSERT(fvec.size() == ideal_point.size());
        GAPP_ASSERT(ideal_point.size() == nadir_point.size());

        FitnessVector fnorm(fvec.size());

        for (size_t i = 0; i < fnorm.size(); i++)
        {
            fnorm[i] = (ideal_point[i] - fvec[i]) / std::max(ideal_point[i] - nadir_point[i], 1E-8);
        }

        return fnorm;
    }


    /* NSGA3 IMPLEMENTATION */

    struct NSGA3::Impl
    {
        struct CandidateTraits
        {
            size_t rank;
            size_t ref_idx;
            double ref_dist;
        };

        RefLineGenerator ref_generator_;
        FitnessMatrix ref_lines_;

        std::vector<CandidateTraits> sol_info_;
        std::vector<size_t> niche_counts_;

        FitnessVector ideal_point_;
        FitnessVector nadir_point_;
        FitnessMatrix extreme_points_;

        /* Generate n reference points in dim dimensions. */
        FitnessMatrix generateReferencePoints(size_t dim, size_t num_points) const;

        /* Update the ideal point approximation using the new points in fmat, assuming maximization. */
        void updateIdealPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Update the extreme points using the new points in fmat, assuming maximization. */
        void updateExtremePoints(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Update the current nadir point based on the extreme points. */
        void updateNadirPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

        /* Recalculate the niche counts of the reference lines based on the ref lines associated with the candidates in pareto_fronts. */
        void recalcNicheCounts(std::span<const FrontElement> pareto_fronts);

        /* Find the closest reference and its distance for each of the points in the fitness matrix. */
        void associatePopWithRefs(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, std::span<const FrontElement> pareto_fronts);

        /* Return true if pop[lhs] is better than pop[rhs]. */
        bool nichedCompare(size_t lhs, size_t rhs) const noexcept;

        /* Return the associated reference direction of a candidate. */
        size_t refIndexOf(const FrontElement& sol) const noexcept;

        /* Return the associated reference direction's distance for a candidate. */
        double refDistOf(const FrontElement& sol) const noexcept;

        /* Return the (unique) reference indices that are associated with at least one element in pareto_fronts, sorted based on their niche counts. */
        std::vector<size_t> referenceSetOf(std::span<const FrontElement> pareto_fronts);

        /* Return the closest solution associated with ref_idx in pareto_fronts. */
        FrontElement& findClosestAssociated(std::span<FrontElement> pareto_fronts, size_t ref_idx) const noexcept;

        /* Increment the niche count of ref, while keeping refs sorted based on niche counts. */
        void incrementNicheCount(std::vector<size_t>& refs, size_t ref);

        /* Create a new population from pareto_fronts. */
        small_vector<size_t> createPopulation(std::span<const FrontElement> pareto_fronts);
    };


    NSGA3::NSGA3(RefLineGenerator gen) :
        pimpl_(std::make_unique<Impl>())
    {
        pimpl_->ref_generator_ = std::move(gen);
    }

    NSGA3::NSGA3(const NSGA3& rhs) :
        pimpl_(std::make_unique<Impl>(*rhs.pimpl_))
    {}

    NSGA3& NSGA3::operator=(const NSGA3& other)
    {
        this->pimpl_ = std::make_unique<Impl>(*other.pimpl_);
        return *this;
    }

    NSGA3& NSGA3::operator=(NSGA3&& other) noexcept
    {
        this->pimpl_ = std::move(other.pimpl_);
        return *this;
    }

    NSGA3::NSGA3(NSGA3&&) noexcept            = default;
    NSGA3::~NSGA3()                           = default;


    FitnessMatrix NSGA3::Impl::generateReferencePoints(size_t dim, size_t num_points) const
    {
        FitnessMatrix ref_points = ref_generator_(dim, num_points);
        std::for_each(ref_points.begin(), ref_points.end(), math::normalizeVector);

        return ref_points;
    }

    void NSGA3::Impl::updateIdealPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);

        detail::elementwise_max(ideal_point_, detail::maxFitness(first, last), detail::inplace_t{});
    }

    void NSGA3::Impl::updateExtremePoints(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);
        GAPP_ASSERT(first->size() == ideal_point_.size());

        const auto popsize = size_t(last - first);

        FitnessMatrix new_extreme_points;
        new_extreme_points.reserve(extreme_points_.nrows(), extreme_points_.ncols());

        for (size_t dim = 0; dim < ideal_point_.size(); dim++)
        {
            const auto weights = weightVector(ideal_point_.size(), dim);
            const auto ASFi = [&](const auto& fvec) { return ASF(ideal_point_, weights, fvec); };
            
            std::vector<double> chebysev_distances(popsize + extreme_points_.size());

            std::transform(first, last, chebysev_distances.begin(), ASFi);
            std::transform(extreme_points_.begin(), extreme_points_.end(), chebysev_distances.begin() + popsize, ASFi);

            const size_t idx = detail::argmin(chebysev_distances.begin(), chebysev_distances.end());

            (idx >= popsize)
                ? new_extreme_points.append_row(extreme_points_[idx - popsize])
                : new_extreme_points.append_row(first[idx]);
        }

        extreme_points_ = std::move(new_extreme_points);
    }

    void NSGA3::Impl::updateNadirPoint(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last)
    {
        updateExtremePoints(first, last);
        nadir_point_ = findFrontNadirPoint(extreme_points_);
    }

    void NSGA3::Impl::recalcNicheCounts(std::span<const FrontElement> pareto_fronts)
    {
        std::fill(niche_counts_.begin(), niche_counts_.end(), 0_sz);
        std::for_each(pareto_fronts.begin(), pareto_fronts.end(), [&](const FrontElement& sol) { niche_counts_[refIndexOf(sol)]++; });
    }

    void NSGA3::Impl::associatePopWithRefs(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last,
                                           std::span<const FrontElement> pareto_fronts)
    {
        GAPP_ASSERT(std::distance(first, last) > 0);
        GAPP_ASSERT(!ref_lines_.empty());

        updateIdealPoint(first, last);
        updateNadirPoint(first, last);

        sol_info_.resize(last - first);

        detail::parallel_for(pareto_fronts.begin(), pareto_fronts.end(), 200, [&](const FrontElement& sol)
        {
            const FitnessVector fnorm = normalizeFitnessVec(first[sol.idx], ideal_point_, nadir_point_);

            auto inverse_distance = [&](const auto& line) { return std::inner_product(fnorm.begin(), fnorm.end(), line.begin(), 0.0); };

            auto closest = detail::max_element(ref_lines_.begin(), ref_lines_.end(), inverse_distance);

            sol_info_[sol.idx].ref_idx  = std::distance(ref_lines_.begin(), closest);
            sol_info_[sol.idx].ref_dist = math::perpendicularDistanceSq(*closest, fnorm);
        });
    }

    bool NSGA3::Impl::nichedCompare(size_t lhs, size_t rhs) const noexcept
    {
        /* This version of the compare implementation is from the U-NSGA-III algorithm. */
        const auto& left  = sol_info_[lhs];
        const auto& right = sol_info_[rhs];

        if (left.ref_idx == right.ref_idx)
        {
            if (left.rank != right.rank)
            {
                return left.rank < right.rank;
            }

            return left.ref_dist < right.ref_dist;
        }

        return rng::randomBool();
    }

    size_t NSGA3::Impl::refIndexOf(const FrontElement& sol) const noexcept
    {
        return sol_info_[sol.idx].ref_idx;
    }

    double NSGA3::Impl::refDistOf(const FrontElement& sol) const noexcept
    {
        return sol_info_[sol.idx].ref_dist;
    }

    std::vector<size_t> NSGA3::Impl::referenceSetOf(std::span<const FrontElement> pareto_fronts)
    {
        std::vector<size_t> refs(pareto_fronts.size());
        std::transform(pareto_fronts.begin(), pareto_fronts.end(), refs.begin(), [this](const FrontElement& sol) { return refIndexOf(sol); });

        detail::erase_duplicates(refs);
        std::sort(refs.begin(), refs.end(), [this](size_t lhs, size_t rhs) { return niche_counts_[lhs] < niche_counts_[rhs]; });

        return refs;
    }

    FrontElement& NSGA3::Impl::findClosestAssociated(std::span<FrontElement> pareto_fronts, size_t ref) const noexcept
    {
        return *detail::min_element(pareto_fronts.begin(), pareto_fronts.end(), [&](const FrontElement& sol)
        {
            return (refIndexOf(sol) == ref) ? refDistOf(sol) : math::inf<double>;
        });
    }

    void NSGA3::Impl::incrementNicheCount(std::vector<size_t>& refs, size_t ref)
    {
        niche_counts_[ref]++;

        auto current = std::find(refs.begin(), refs.end(), ref);
        auto first_eq = std::find_if(std::next(current), refs.end(), [&](size_t idx) { return niche_counts_[idx] >= niche_counts_[ref]; });

        std::iter_swap(current, std::prev(first_eq));
    }

    small_vector<size_t> NSGA3::Impl::createPopulation(std::span<const FrontElement> pareto_fronts)
    {
        small_vector<size_t> new_pop;
        std::vector<Impl::CandidateTraits> new_traits;

        new_pop.reserve(pareto_fronts.size());
        new_traits.reserve(pareto_fronts.size());

        for (const FrontElement& sol : pareto_fronts)
        {
            new_pop.push_back(sol.idx);
            new_traits.push_back(sol_info_[sol.idx]);
        }

        sol_info_ = std::move(new_traits);

        return new_pop;
    }


    /* NSGA3 INTERFACE FUNCTIONS */

    void NSGA3::initializeImpl(const GaInfo& ga)
    {
        GAPP_ASSERT(ga.population_size() != 0);
        GAPP_ASSERT(ga.num_objectives() > 1, "The number of objectives must be greater than 1 for the NSGA-III algorithm.");

        const FitnessMatrix& fitness_matrix = ga.fitness_matrix();

        pimpl_->ideal_point_ = detail::maxFitness(fitness_matrix.begin(), fitness_matrix.end());
        pimpl_->extreme_points_ = {};

        pimpl_->ref_lines_ = pimpl_->generateReferencePoints(ga.num_objectives(), ga.population_size());
        pimpl_->niche_counts_.resize(pimpl_->ref_lines_.size());

        ParetoFronts pareto_fronts = nonDominatedSort(fitness_matrix.begin(), fitness_matrix.end());

        pimpl_->sol_info_.resize(ga.population_size());
        std::for_each(pareto_fronts.begin(), pareto_fronts.end(), [this](const FrontElement& sol) { pimpl_->sol_info_[sol.idx].rank = sol.rank; });

        pimpl_->associatePopWithRefs(fitness_matrix.begin(), fitness_matrix.end(), pareto_fronts);
        pimpl_->recalcNicheCounts(pareto_fronts);
    }

    small_vector<size_t> NSGA3::nextPopulationImpl(const GaInfo& ga, const PopulationView& pop)
    {
        GAPP_ASSERT(ga.num_objectives() > 1);

        const size_t popsize = ga.population_size();
        const FitnessMatrix fitness_matrix = detail::toFitnessMatrix(pop);

        auto pareto_fronts = nonDominatedSort(fitness_matrix);
        auto partial_front = pareto_fronts.partialFront(popsize);

        pimpl_->sol_info_.resize(pop.size());
        std::for_each(pareto_fronts.begin(), pareto_fronts.end(), [this](const FrontElement& sol) { pimpl_->sol_info_[sol.idx].rank = sol.rank; });

        /* The ref lines of the candidates after partial_last are irrelevant, as they can never be part of the next population. */
        pimpl_->associatePopWithRefs(fitness_matrix.begin(), fitness_matrix.end(), { pareto_fronts.begin(), partial_front.end() });
        /* The niche counts should be calculated excluding the partial front for now. */
        pimpl_->recalcNicheCounts({ pareto_fronts.begin(), partial_front.begin() });

        /* Find the reference lines associated with the partial front. */
        std::vector<size_t> ref_indices = pimpl_->referenceSetOf(partial_front);

        /* Move the best elements to the front of the partial front. */
        for (size_t i = 0; i < popsize - (partial_front.begin() - pareto_fronts.begin()); i++)
        {
            const size_t min_niche_count = pimpl_->niche_counts_[ref_indices.front()];
            const auto last_minimal_ref = std::find_if(ref_indices.begin(), ref_indices.end(), [&](size_t idx)
            {
                return pimpl_->niche_counts_[idx] != min_niche_count;
            });

            const size_t selected_ref_idx = *rng::randomElement(ref_indices.begin(), last_minimal_ref);
            pimpl_->incrementNicheCount(ref_indices, selected_ref_idx);

            const size_t associated_sol_count = std::count_if(pareto_fronts.begin(), pareto_fronts.end(), [&](const FrontElement& sol)
            {
                return pimpl_->refIndexOf(sol) == selected_ref_idx;
            });

            FrontElement& closest = pimpl_->findClosestAssociated(partial_front, selected_ref_idx);

            /* Move the selected candidate to the front of the partial front so it can't be selected again. */
            std::swap(closest, partial_front[i]);

            /* If the selected candidate was the only one in the partial front associated with this reference direction,
               the reference direction needs to also be removed. */
            if (associated_sol_count == 1) detail::erase_first_stable(ref_indices, selected_ref_idx);
        }

        pareto_fronts.resize(popsize);

        return pimpl_->createPopulation(pareto_fronts);
    }

    size_t NSGA3::selectImpl(const GaInfo&, const PopulationView& pop) const
    {
        GAPP_ASSERT(!pop.empty());

        const size_t idx1 = rng::randomIndex(pop);
        const size_t idx2 = rng::randomIndex(pop);

        return pimpl_->nichedCompare(idx1, idx2) ? idx1 : idx2;
    }

    small_vector<size_t> NSGA3::optimalSolutionsImpl(const GaInfo&, const PopulationView&) const
    {
        return detail::find_indices(pimpl_->sol_info_, [](const auto& sol) { return sol.rank == 0; });
    }

} // namespace gapp::algorithm
