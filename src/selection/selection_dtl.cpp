/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "selection_dtl.hpp"
#include "../utility/rng.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include "../utility/algorithm.hpp"
#include "../utility/functional.hpp"
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>

namespace genetic_algorithm::selection::dtl
{
    std::vector<double> rouletteWeights(const FitnessMatrix& fmat)
    {
        auto fvec = detail::toFitnessVector(fmat);

        /* Roulette selection wouldn't work for negative fitness values. */
        double offset = *std::min_element(fvec.begin(), fvec.end());
        offset = 2.0 * offset;              /* The selection probability of the worst candidate should also be > 0. */
        offset = std::min(0.0, offset);     /* Only adjust fitness values if it's neccesary (there are negative fitness values). */

        std::transform(fvec.begin(), fvec.end(), fvec.begin(),
        [offset](double f)
        {
            return f - offset;
        });

        return fvec;
    }

    std::vector<double> rankWeights(const FitnessMatrix& fmat, double wmin, double wmax)
    {
        assert(0.0 <= wmin && wmin <= wmax);

        auto fvec = detail::toFitnessVector(fmat);
        auto indices = detail::argsort(fvec.begin(), fvec.end());

        std::vector<double> weights(fmat.size());
        for (size_t i = 0; i < indices.size(); i++)
        {
            double t = i / (weights.size() - 1.0);
            weights[indices[i]] = std::lerp(wmin, wmax, t);
        }

        return weights;
    }

    std::vector<double> sigmaWeights(const FitnessMatrix& fmat, double scale)
    {
        assert(scale > 1.0);

        auto fvec = detail::toFitnessVector(fmat);
        double fmean = detail::mean(fvec);
        double fdev = std::max(detail::stdDev(fvec, fmean), 1E-6);

        std::transform(fvec.begin(), fvec.end(), fvec.begin(),
        [fmean, fdev, scale](double f)
        {
            double weight = 1.0 + (f - fmean) / (scale * fdev);

            return std::max({ weight, 0.0 });  /* If ( fitness < (f_mean - scale * SD) ) the weight could be negative. */
        });

        return fvec;
    }

    std::vector<double> boltzmannWeights(const FitnessMatrix& fmat, double temperature)
    {
        auto fvec = detail::toFitnessVector(fmat);
        auto [fmin, fmax] = std::minmax_element(fvec.begin(), fvec.end());

        std::transform(fvec.begin(), fvec.end(), fvec.begin(),
        // dont try to capture the iterators by ref or value here
        [fmin = *fmin, fmax = *fmax, temperature](double f) noexcept
        {
            double df = std::max(fmax - fmin, 1E-6);
            double fnorm = (f - fmin) / df;

            return std::exp(fnorm / temperature);
        });

        return fvec;
    }

    double boltzmannDefaultTemp(size_t gen, size_t max_gen) noexcept
    {
        return -4.0 / (1.0 + std::exp(-10.0 * (double(gen) / max_gen) + 3.0)) + 4.0 + 0.25;
    }

    std::vector<double> weightsToCdf(const std::vector<double>& weights)
    {
        double wmean = detail::mean(weights);

        return detail::map(weights,
        [cdf = 0.0, wmean, n = weights.size()](double w) mutable
        {
            return cdf += w / wmean / n;
        });
    }

    ParetoFronts nonDominatedSort(const FitnessMatrix& fmat)
    {
        size_t pop_size = fmat.size();

        /* Find the number of candidates which dominate each candidate, and the indices of the candidates it dominates. */
        std::vector<size_t> better_count(pop_size, 0);
        std::vector<std::vector<size_t>> worse_indices(pop_size);

        std::for_each(worse_indices.begin(), worse_indices.end(), [&pop_size](std::vector<size_t>& vec) { vec.reserve(pop_size); });

        for (size_t lhs = 0; lhs < pop_size; lhs++)
        {
            for (size_t rhs = 1; rhs < lhs; rhs++)
            {
                if (detail::paretoCompareLess(fmat[lhs], fmat[rhs]))
                {
                    better_count[lhs]++;
                    worse_indices[rhs].push_back(lhs);
                }
                else if (detail::paretoCompareLess(fmat[rhs], fmat[lhs])) [[likely]]
                {
                    better_count[rhs]++;
                    worse_indices[lhs].push_back(rhs);
                }
            }
        }

        /* [idx, rank] */
        ParetoFronts sorted_indices;
        sorted_indices.reserve(pop_size);

        /* Find the indices of all non-dominated candidates (the first/best pareto front). */
        for (size_t i = 0; i < pop_size; i++)
        {
            if (better_count[i] == 0)
            {
                sorted_indices.emplace_back(i, 0);
            }
        }

        /* Find all the other pareto fronts. */
        auto current_first = sorted_indices.cbegin();
        auto current_last  = sorted_indices.cend();

        while (sorted_indices.size() != pop_size)
        {
            size_t next_rank = current_first->second + 1;

            /* Remove the current front from the population and find the next one. */
            for (; current_first != current_last; current_first++)
            {
                for (const auto& worse_idx : worse_indices[current_first->first])
                {
                    if (--better_count[worse_idx] == 0)
                    {
                        current_last--;
                        sorted_indices.emplace_back(worse_idx, next_rank);
                        current_last++;
                    }
                }
            }
            current_last = sorted_indices.cend();
        }

        return sorted_indices;
    }

    std::vector<size_t> paretoRanks(const ParetoFronts& pareto_fronts)
    {
        std::vector<size_t> ranks(pareto_fronts.size());

        for (const auto& [idx, rank] : pareto_fronts)
        {
            ranks[idx] = rank;
        }

        return ranks;
    }

    ParetoFronts::iterator nextFrontBegin(ParetoFronts::iterator current, ParetoFronts::iterator last) noexcept
    {
        return std::find_if(current, last,
        [current_rank = current->second](const std::pair<size_t, size_t>& elem) noexcept
        {
            size_t rank = elem.second;
            return rank != current_rank;
        });
    }

    auto paretoFrontBounds(ParetoFronts& pareto_fronts) -> std::vector<std::pair<ParetoFronts::iterator, ParetoFronts::iterator>>
    {
        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        std::vector<IterPair> front_bounds;
        front_bounds.reserve(pareto_fronts.back().second);

        for (auto first = pareto_fronts.begin(); first != pareto_fronts.end(); )
        {
            auto last = nextFrontBegin(first, pareto_fronts.end());
            front_bounds.emplace_back(first, last);
            first = last;
        }

        return front_bounds;
    }

    std::vector<double> crowdingDistances(const FitnessMatrix& fmat, ParetoFronts pfronts)
    {
        std::vector<double> distances(fmat.size(), 0.0);

        using Iter = ParetoFronts::iterator;
        using IterPair = std::pair<Iter, Iter>;

        auto front_bounds = paretoFrontBounds(pfronts);

        std::for_each(front_bounds.begin(), front_bounds.end(),
        [&distances, &fmat](const IterPair& bounds)
        {
            const auto& [first, last] = bounds;

            /* Calculate distances along each fitness dimension */
            for (size_t dim = 0; dim < fmat[0].size(); dim++)
            {
                std::sort(first, last,
                [&fmat, dim](const auto& lhs, const auto& rhs) noexcept
                {
                    return fmat[lhs.first][dim] < fmat[rhs.first][dim];
                });

                const auto& front = first->first;
                const auto& back  = (last - 1)->first;

                double finterval = fmat[back][dim] - fmat[front][dim];
                finterval = std::max(finterval, 1E-6);

                distances[front] = std::numeric_limits<double>::infinity();
                distances[back]  = std::numeric_limits<double>::infinity();

                for (auto it = std::next(first); it < std::prev(last); it++)
                {
                    size_t this_ = it->first;
                    size_t next = std::next(it)->first;
                    size_t prev = std::prev(it)->first;

                    distances[this_] += (fmat[next][dim] - fmat[prev][dim]) / finterval;
                }
            }
        });

        return distances;
    }

    std::vector<Point> generateRefPoints(size_t n, size_t dim)
    {
        using namespace std;
        assert(n > 0);
        assert(dim > 1);

        /* Generate reference point candidates randomly. */
        size_t ratio = max(size_t{ 10 }, 2 * dim);
        vector<Point> candidates(ratio * n - 1);
        generate(candidates.begin(), candidates.end(), [&dim]() { return rng::randomSimplexPoint(dim); });

        vector<Point> refs;
        refs.reserve(n);
        refs.push_back(rng::randomSimplexPoint(dim));

        vector<double> min_distances(candidates.size(), numeric_limits<double>::infinity());
        while (refs.size() < n)
        {
            /* Calc the distance of each candidate to the closest ref point. */
            transform(GA_EXECUTION_UNSEQ, candidates.begin(), candidates.end(), min_distances.begin(), min_distances.begin(),
            [&refs](const Point& candidate, double dmin) noexcept
            {
                double d = detail::euclideanDistanceSq(candidate, refs.back());
                return min(dmin, d);
            });

            /* Add the candidate with highest min_distance to the refs. */
            size_t argmax = detail::argmax(min_distances.begin(), min_distances.end());
            refs.push_back(move(candidates[argmax]));

            /* Remove the added candidate and the corresponding min_distance. */
            swap(candidates[argmax], candidates.back());
            candidates.pop_back();
            swap(min_distances[argmax], min_distances.back());
            min_distances.pop_back();
        }

        return refs;
    }

    std::pair<size_t, double> findClosestRef(const std::vector<RefPoint>& refs, const Point& p)
    {
        assert(!refs.empty());

        auto distances = detail::map(refs,
        [&p](const RefPoint& ref) noexcept
        {
            return detail::perpendicularDistanceSq(ref.point, p);
        });

        auto idx = detail::argmin(distances.begin(), distances.end());

        return { idx, distances[idx] };
    }

    auto ASF(std::vector<double> z, std::vector<double> w) noexcept
        -> std::function<double(const std::vector<double>&)>
    {
        assert(!w.empty());
        assert(w.size() == z.size());

        return [z = std::move(z), w = std::move(w)](const std::vector<double>& f) noexcept
        {
            double dmax = std::abs(f[0] - z[0]) / w[0];
            for (size_t i = 0; i < f.size(); i++)
            {
                dmax = std::max(dmax, std::abs(f[i] - z[i]) / w[i]);
            }
            return dmax;
        };
    }

} // namespace genetic_algorithm::selection::dtl