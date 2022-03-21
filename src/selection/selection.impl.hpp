/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_IMPL_HPP
#define GA_SELECTION_IMPL_HPP

#include "../population.hpp"
#include "selection_dtl.hpp"
#include "../rng.hpp"
#include "../algorithms/ga_base.decl.hpp"

#include <vector>
#include <algorithm>
#include <iterator>
#include <limits>
#include <exception>
#include <cassert>

namespace genetic_algorithm::selection
{
    template<gene T>
    inline void Roulette<T>::prepare(const GA<T>&, const Population<T>& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate<T>& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        auto selection_weights = dtl::rouletteWeights(fitnessVector(pop));

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    template<gene T>
    inline Candidate<T> Roulette<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(!pop.empty() && pop.size() == cdf_.size());

        return pop[dtl::sampleCdf(cdf_)];
    }

    template<gene T>
    inline Tournament<T>::Tournament(size_t size)
    {
        this->size(size);
    }

    template<gene T>
    inline Tournament<T>::Tournament(const GA<T>&, size_t size)
    {
        this->size(size);
    }

    template<gene T>
    inline void Tournament<T>::size(size_t size)
    {
        if (size < 2) { throw std::invalid_argument("The tournament size must be at least 2."); }

        tourney_size_ = size;
    }

    template<gene T>
    inline void Tournament<T>::prepare(const GA<T>&, const Population<T>&)
    { /* Nothing to do for tournament selection. */ }

    template<gene T>
    inline Candidate<T> Tournament<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(pop.size() >= tourney_size_);
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate<T>& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        auto candidates = rng::sampleUnique(pop.size(), tourney_size_);

        size_t idx = *std::max_element(candidates.begin(), candidates.end(),
        [&pop](size_t lidx, size_t ridx)
        {
            return pop[lidx].fitness[0] < pop[ridx].fitness[0];
        });

        return pop[idx];
    }

    template<gene T>
    inline Rank<T>::Rank(double min_weight, double max_weight)
    {
        this->weights(min_weight, max_weight);
    }

    template<gene T>
    inline Rank<T>::Rank(const GA<T>&, double min_weight, double max_weight)
    {
        this->weights(min_weight, max_weight);
    }

    template<gene T>
    inline void Rank<T>::min_weight(double min_weight)
    {
        this->weights(min_weight, max_weight_);
    }

    template<gene T>
    inline void Rank<T>::max_weight(double max_weight)
    {
        this->weights(min_weight_, max_weight);
    }

    template<gene T>
    inline void Rank<T>::weights(double min_weight, double max_weight)
    {
        if (!(0.0 <= min_weight && min_weight <= max_weight))
        {
            throw std::invalid_argument("The minimum weight must be in the closed interval [0.0, max_weight].");
        }
        if (!(min_weight <= max_weight && max_weight <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("The maximum weight must be in the closed interval [min_weight, DBL_MAX].");
        }

        min_weight_ = min_weight;
        max_weight_ = max_weight;
    }

    template<gene T>
    inline void Rank<T>::prepare(const GA<T>&, const Population<T>& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate<T>& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        auto selection_weights = dtl::rankWeights(fitnessVector(pop), min_weight_, max_weight_);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    template<gene T>
    inline Candidate<T> Rank<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(!pop.empty() && pop.size() == cdf_.size());

        return pop[dtl::sampleCdf(cdf_)];
    }

    template<gene T>
    inline Sigma<T>::Sigma(double scale)
    {
        this->scale(scale);
    }

    template<gene T>
    inline Sigma<T>::Sigma(const GA<T>&, double scale)
    {
        this->scale(scale);
    }

    template<gene T>
    inline void Sigma<T>::scale(double scale)
    {
        if (!(1.0 <= scale && scale <= std::numeric_limits<double>::max()))
        {
            throw std::invalid_argument("Scale must be in the closed interval [1.0, DBL_MAX].");
        }

        scale_ = scale;
    }

    template<gene T>
    inline void Sigma<T>::prepare(const GA<T>&, const Population<T>& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate<T>& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        auto selection_weights = dtl::sigmaWeights(fitnessVector(pop), scale_);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    template<gene T>
    inline Candidate<T> Sigma<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(pop.size() == cdf_.size());

        return pop[dtl::sampleCdf(cdf_)];
    }

    template<gene T>
    inline Boltzmann<T>::Boltzmann(TemperatureFunction f) : temperature_(std::move(f))
    {}

    template<gene T>
    inline Boltzmann<T>::Boltzmann(const GA<T>&, TemperatureFunction f) : temperature_(std::move(f))
    {}



    template<gene T>
    inline void Boltzmann<T>::prepare(const GA<T>& ga, const Population<T>& pop)
    {
        assert(std::all_of(pop.begin(), pop.end(), [](const Candidate<T>& sol) { return sol.fitness.size() == 1 && sol.is_evaluated; }));

        double temp = std::invoke(temperature_, ga.generation_cntr(), ga.max_gen());
        auto selection_weights = dtl::boltzmannWeights(fitnessVector(pop), temp);

        cdf_ = dtl::weightsToCdf(selection_weights);
    }

    template<gene T>
    inline Candidate<T> Boltzmann<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(!pop.empty() && pop.size() == cdf_.size());

        return pop[dtl::sampleCdf(cdf_)];
    }

    template<gene T>
    inline void NSGA2<T>::init(const GA<T>& ga)
    {
        assert(ga.num_objectives() > 1);
        assert(!ga.population().empty());

        auto pfronts = dtl::nonDominatedSort(fitnessMatrix(ga.population()));
        ranks_ = pfronts.ranks;
        dists_ = dtl::crowdingDistances(fitnessMatrix(ga.population()), pfronts.idxs);
    }

    template<gene T>
    inline void NSGA2<T>::prepare(const GA<T>&, const Population<T>&)
    { /* Nothing to do, the ranks and distances from the previous nextPopulation call are fine */ }

    template<gene T>
    inline Candidate<T> NSGA2<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop.size());
        size_t idx2 = rng::randomIdx(pop.size());

        return crowdedCompare(idx1, idx2) ? pop[idx1] : pop[idx2];
    }

    template<gene T>
    inline Population<T> NSGA2<T>::nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children)
    {
        assert(!children.empty());

        Population<T> new_pop;
        new_pop.reserve(ga.population_size());

        old_pop.insert(old_pop.end(), std::make_move_iterator(children.begin()),
                                      std::make_move_iterator(children.end()));

        assert(std::all_of(old_pop.begin(), old_pop.end(), [](const Candidate<T>& sol) { return sol.is_evaluated; }));

        auto pfronts = dtl::nonDominatedSort(fitnessMatrix(old_pop));
        ranks_ = pfronts.ranks;
        dists_ = dtl::crowdingDistances(fitnessMatrix(old_pop), pfronts.idxs);

        /* Keep track of the ranks of the candidates added to new_pop to avoid a second sort. */
        std::vector<size_t> new_ranks;
        new_ranks.reserve(ga.population_size());
        std::vector<double> new_dists;
        new_dists.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        while (new_pop.size() + pfronts.idxs[front_idx].size() <= ga.population_size())
        {
            for (const auto& idx : pfronts.idxs[front_idx])
            {
                new_pop.push_back(std::move(old_pop[idx]));
                new_ranks.push_back(ranks_[idx]);
                new_dists.push_back(dists_[idx]);
            }
            front_idx++;
        }

        /* Add the remaining candidates from the partial front if there is one. */
        if (new_pop.size() != ga.population_size())
        {
            /* Keep track of the candidates added from the partial front, 
             * the crowding distances will only need to be updated for these candidates. */
            std::vector<size_t> added_indices(ga.population_size() - new_pop.size());
            std::iota(added_indices.begin(), added_indices.end(), new_pop.size());

            std::vector<size_t> partial_front = pfronts.idxs[front_idx];

            std::sort(partial_front.begin(), partial_front.end(),
            [this](size_t lidx, size_t ridx)
            {
                return crowdedCompare(lidx, ridx);
            });

            for (const auto& idx : partial_front)
            {
                new_pop.push_back(std::move(old_pop[idx]));
                new_ranks.push_back(ranks_[idx]);
                new_dists.push_back(dists_[idx]);

                if (new_pop.size() == ga.population_size()) break;
            }

            auto changed_dists = dtl::crowdingDistances(fitnessMatrix(new_pop), { added_indices });
            for (size_t idx : added_indices)
            {
                new_dists[idx] = changed_dists[idx];
            }
            
        }
        ranks_ = std::move(new_ranks);
        dists_ = std::move(new_dists);

        return new_pop;
    }

    template<gene T>
    inline bool NSGA2<T>::crowdedCompare(size_t lidx, size_t ridx) const
    {
        if (ranks_[lidx] != ranks_[ridx])
        {
            return ranks_[lidx] < ranks_[ridx];
        }
        else
        {
            return dists_[lidx] > dists_[ridx];
        }
    }

    template<gene T>
    inline void NSGA3<T>::init(const GA<T>& ga)
    {
        assert(ga.num_objectives() > 1);
        assert(!ga.population().empty());

        ref_points_ = dtl::generateRefPoints(ga.population_size(), ga.num_objectives());

        ideal_point_ = std::vector(ga.num_objectives(), -std::numeric_limits<double>::max());
        updateIdealPoint(ideal_point_, ga.population());
        extreme_points_ = initExtremePoints(ga.population(), ideal_point_);
        nadir_point_ = nadirPoint(extreme_points_);

        sol_props_ = std::vector<CandidateInfo>(ga.population_size());
        associatePopWithRefs(ga.population());
        ref_niche_counts_ = calcNicheCounts(ga, ga.population(), sol_props_);

        auto pfronts = dtl::nonDominatedSort(fitnessMatrix(ga.population()));
        for (size_t i = 0; i < sol_props_.size(); i++)
        {
            sol_props_[i].rank = pfronts.ranks[i];
        }
    }

    template<gene T>
    inline void NSGA3<T>::updateIdealPoint(Point& ideal_point, const Population<T>& pop)
    {
        for (const auto& sol : pop)
        {
            for (size_t i = 0; i < ideal_point.size(); i++)
            {
                ideal_point[i] = std::max(ideal_point[i], sol.fitness[i]);
            }
        }
    }

    template<gene T>
    inline auto NSGA3<T>::initExtremePoints(const Population<T>& pop, const Point& ideal_point) -> std::vector<Point>
    {
        assert(!pop.empty());

        std::vector<Point> extreme_points(ideal_point.size(), Point(ideal_point.size()));

        /* Identify and update extreme points for each objective axis. */
        for (size_t i = 0; i < extreme_points.size(); i++)
        {
            std::vector<double> weights(extreme_points[i].size(), 1E-6);
            weights[i] = 1.0;

            /* Find the solution or extreme point with the lowest Chebysev distance to the objective axis. */
            double dmin = std::numeric_limits<double>::infinity();
            std::vector<double> extreme_point = pop[0].fitness;
            for (const auto& sol : pop)
            {
                double dist = dtl::ASF(sol.fitness, ideal_point, weights);

                if (dist < dmin)
                {
                    dmin = dist;
                    extreme_point = sol.fitness;
                }
            }

            extreme_points[i] = extreme_point;
        }

        return extreme_points;
    }

    template<gene T>
    inline void NSGA3<T>::updateExtremePoints(std::vector<Point>& extreme_points, const Population<T>& pop, const Point& ideal_point)
    {
        assert(!pop.empty());

        /* Identify and update extreme points for each objective axis. */
        for (size_t i = 0; i < extreme_points.size(); i++)
        {
            std::vector<double> weights(extreme_points[i].size(), 1E-6);
            weights[i] = 1.0;

            /* Find the solution or extreme point with the lowest Chebysev distance to the objective axis. */
            double dmin = std::numeric_limits<double>::infinity();
            std::vector<double> new_extreme_point = extreme_points[i];

            for (const auto& sol : pop)
            {
                double dist = dtl::ASF(sol.fitness, ideal_point, weights);

                if (dist < dmin)
                {
                    dmin = dist;
                    new_extreme_point = sol.fitness;
                }
            }
            for (const auto& old_extreme_point : extreme_points)
            {
                double dist = dtl::ASF(old_extreme_point, ideal_point, weights);

                if (dist < dmin)
                {
                    dmin = dist;
                    new_extreme_point = old_extreme_point;
                }
            }

            extreme_points[i] = new_extreme_point;
        }
    }

    template<gene T>
    inline auto NSGA3<T>::nadirPoint(const std::vector<Point>& extreme_points) -> Point
    {
        Point nadir(extreme_points.size());

        /* Nadir point estimate = minimum of extreme points along each objective axis. */
        for (size_t i = 0; i < extreme_points.size(); i++)
        {
            nadir[i] = extreme_points[0][i];
            for (size_t j = 1; j < extreme_points.size(); j++)
            {
                nadir[i] = std::min(nadir[i], extreme_points[j][i]);
            }
        }

        return nadir;
    }

    template<gene T>
    inline void NSGA3<T>::prepare(const GA<T>&, const Population<T>&)
    { /* Nothing to do */ }

    template<gene T>
    inline bool NSGA3<T>::nichedCompare(size_t lidx, size_t ridx) const
    {
        if (sol_props_[lidx].rank != sol_props_[ridx].rank)
        {
            return sol_props_[lidx].rank < sol_props_[ridx].rank;
        }
        else if (sol_props_[lidx].niche_count != sol_props_[ridx].niche_count)
        {
            return sol_props_[lidx].niche_count < sol_props_[ridx].niche_count;
        }
        else
        {
            return sol_props_[lidx].ref_dist < sol_props_[ridx].ref_dist;
        }
    }

    template<gene T>
    inline Candidate<T> NSGA3<T>::select(const GA<T>&, const Population<T>& pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop.size());
        size_t idx2 = rng::randomIdx(pop.size());

        return nichedCompare(idx1, idx2) ? pop[idx1] : pop[idx2];
    }

    template<gene T>
    inline void NSGA3<T>::associatePopWithRefs(const Population<T>& pop)
    {
        assert(!pop.empty());
        assert(std::all_of(pop.begin(), pop.end(), [&pop](const Candidate<T>& sol) { return sol.fitness.size() == pop[0].fitness.size(); }));

        auto fnorms = fitnessMatrix(pop);

        std::for_each(GA_EXECUTION_UNSEQ, fnorms.begin(), fnorms.end(),
        [this](std::vector<double>& fnorm)
        {
            for (size_t i = 0; i < fnorm.size(); i++)
            {
                fnorm[i] -= ideal_point_[i];
                fnorm[i] /= std::min(nadir_point_[i] - ideal_point_[i], -1E-6);
            }
        });

        /* Associate each candidate with the closest reference point. */
        sol_props_.resize(pop.size());
        std::transform(GA_EXECUTION_UNSEQ, fnorms.begin(), fnorms.end(), sol_props_.begin(), sol_props_.begin(),
        [this](const std::vector<double>& f, CandidateInfo& info) -> CandidateInfo
        {
            std::tie(info.ref_idx, info.ref_dist) = dtl::findClosestRef(ref_points_, f);

            return info;
        });
    }

    template<gene T>
    inline std::vector<size_t> NSGA3<T>::calcNicheCounts(const GA<T>& ga, const Population<T>& pop, std::vector<CandidateInfo>& props)
    {
        assert(pop.size() == props.size());

        /* Calculate niche counts for each of the reference points */
        std::vector<size_t> ref_niche_counts(ga.population_size(), 0U);
        for (const auto& info : props)
        {
            ref_niche_counts[info.ref_idx]++;
        }

        /* Assign the niche counts to the candidates too. */
        for (auto& info : props)
        {
            info.niche_count = ref_niche_counts[info.ref_idx];
        }

        return ref_niche_counts;
    }

    template<gene T>
    inline Population<T> NSGA3<T>::nextPopulation(const GA<T>& ga, Population<T>& old_pop, CandidateVec<T>& children)
    {
        assert(!children.empty());

        updateIdealPoint(ideal_point_, old_pop);
        updateExtremePoints(extreme_points_, old_pop, ideal_point_);
        nadir_point_ = nadirPoint(extreme_points_);

        Population<T> new_pop;
        new_pop.reserve(ga.population_size());

        old_pop.insert(old_pop.end(), std::make_move_iterator(children.begin()),
                                      std::make_move_iterator(children.end()));
        sol_props_.resize(old_pop.size());

        assert(std::all_of(old_pop.begin(), old_pop.end(), [](const Candidate<T>& sol) { return sol.is_evaluated; }));

        auto pfronts = dtl::nonDominatedSort(fitnessMatrix(old_pop));
        for (size_t i = 0; i < sol_props_.size(); i++)
        {
            sol_props_[i].rank = pfronts.ranks[i];
        }
        associatePopWithRefs(old_pop);

        std::vector<CandidateInfo> new_props;
        new_props.reserve(ga.population_size());

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        while (new_pop.size() + pfronts.idxs[front_idx].size() <= ga.population_size())
        {
            for (const auto& idx : pfronts.idxs[front_idx])
            {
                new_pop.push_back(std::move(old_pop[idx]));
                new_props.push_back(sol_props_[idx]);
            }
            front_idx++;
        }
        std::vector<size_t> ref_niche_counts = calcNicheCounts(ga, new_pop, new_props);

        /* Add remaining candidates from the partial front if there is one. */
        std::vector<size_t> partial_front = pfronts.idxs[front_idx];
        while (new_pop.size() != ga.population_size())
        {
            /* Find the lowest niche count in the partial front. */
            size_t min_count = std::numeric_limits<size_t>::max();
            for (const auto& idx : partial_front)
            {
                min_count = std::min(min_count, ref_niche_counts[sol_props_[idx].ref_idx]);
            }

            /* Find the reference points with minimal niche counts, and pick one. */
            std::vector<size_t> refs;
            for (const auto& idx : partial_front)
            {
                size_t ref = sol_props_[idx].ref_idx;
                if (ref_niche_counts[ref] == min_count && std::find(refs.begin(), refs.end(), ref) == refs.end())
                {
                    refs.push_back(ref);
                }
            }
            size_t ref = refs[rng::randomIdx(refs.size())];

            /* Find the idx of the closest sol in the partial front associated with this ref point. */
            size_t sol_idx = partial_front[0];
            double min_distance = std::numeric_limits<double>::max();
            for (const auto& idx : partial_front)
            {
                if (sol_props_[idx].ref_idx == ref && sol_props_[idx].ref_dist < min_distance)
                {
                    min_distance = sol_props_[idx].ref_dist;
                    sol_idx = idx;
                }
            }

            /* Move this candidate to new_pop and increment the associated niche count. */
            new_pop.push_back(std::move(old_pop[sol_idx]));
            new_props.push_back(sol_props_[sol_idx]);
            partial_front.erase(std::remove(partial_front.begin(), partial_front.end(), sol_idx), partial_front.end());

            ref_niche_counts[ref]++;
        }

        ref_niche_counts_ = calcNicheCounts(ga, new_pop, new_props);
        sol_props_ = new_props;

        return new_pop;
    }

} // namespace genetic_algorithm::selection

#endif // !GA_SELECTION_IMPL_HPP