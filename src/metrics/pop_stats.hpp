/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_METRICS_POP_STATS_HPP
#define GA_METRICS_POP_STATS_HPP

#include "../population/candidate.hpp"
#include <span>

namespace genetic_algorithm::detail
{
    /* Return the minimum fitness values of a fitness matrix along each objective axis. */
    FitnessVector minFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the maximum fitness values of a fitness matrix along each objective axis. */
    FitnessVector maxFitness(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the mean fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessMean(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the variance of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the variance of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessVariance(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, std::span<const double> fitness_mean);

    /* Return the standard deviation of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last);

    /* Return the standard deviation of the fitness values of a fitness matrix along each objective axis. */
    FitnessVector fitnessStdDev(FitnessMatrix::const_iterator first, FitnessMatrix::const_iterator last, std::span<const double> fitness_mean);

    /* Compute the hypervolume of a set of points relative to a reference point. Works for any number of dimensions. */
    double hypervolume(const FitnessMatrix& fmat, std::span<const double> ref_point);

} // namespace genetic_algorithm::detail

#endif // !GA_METRICS_POP_STATS_HPP
