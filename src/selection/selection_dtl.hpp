/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_SELECTION_DTL_HPP
#define GA_SELECTION_DTL_HPP

#include "../population/population.hpp"
#include <vector>
#include <utility>
#include <cstddef>

namespace genetic_algorithm::selection::dtl
{
    using genetic_algorithm::detail::FitnessVector;
    using genetic_algorithm::detail::FitnessMatrix;

    /* Calculate the selection weights of the population for the roulette selection. */
    std::vector<double> rouletteWeights(const FitnessMatrix& pop);

    /* Calculate the selection weights of the population for the rank selection. */
    std::vector<double> rankWeights(const FitnessMatrix& pop, double wmin, double wmax);

    /* Calculate the selection weights of the population for the sigma selection. */
    std::vector<double> sigmaWeights(const FitnessMatrix& pop, double scale);

    /* Calculate the selection weights of the population for the Boltzmann selection. */
    std::vector<double> boltzmannWeights(const FitnessMatrix& pop, double temperature);

    /* Default temperature function used for the Boltzmann selection. */
    double boltzmannDefaultTemp(size_t gen, size_t max_gen);

    /* Calculate the cumulative distribution function of the population from the selection weights. */
    std::vector<double> weightsToCdf(const std::vector<double>& selection_weights);


    struct ParetoFrontsInfo
    {
    public:
        std::vector<std::vector<size_t>> idxs;
        std::vector<size_t> ranks;

        ParetoFrontsInfo(std::vector<std::vector<size_t>> idxs, std::vector<size_t> ranks)
            : idxs(std::move(idxs)), ranks(std::move(ranks)) {}
    };

    /* Non-dominated sorting for the multi-objective algorithms. Returns the pareto fronts of the population and the ranks of each candidate. */
    ParetoFrontsInfo nonDominatedSort(const FitnessMatrix& fmat);

    /* Calculate the crowding distances of the solutions in the NSGA2 algorithm. */
    std::vector<double> crowdingDistances(const FitnessMatrix& fmat, std::vector<std::vector<size_t>>& pfronts);

    /* Calculate the crowding distances of the solutions in the NSGA2 algorithm. */
    std::vector<double> crowdingDistances(const FitnessMatrix& fmat, std::vector<std::vector<size_t>>&& pfronts);

    /* Generate n reference points on the unit simplex in dim dimensions (for the NSGA-III algorithm). */
    std::vector<std::vector<double>> generateRefPoints(size_t n, size_t dim);

    /* Find the index and distance of the closest reference line to the point p. */
    std::pair<size_t, double> findClosestRef(const std::vector<std::vector<double>>& refs, const std::vector<double>& p);

    /* Achievement scalarization function for the NSGA-III algorithm. */
    double ASF(const std::vector<double>& f, const std::vector<double>& z, const std::vector<double>& w);

} // namespace genetic_algorithm::selection::dtl

#endif // !GA_SELECTION_DTL_HPP