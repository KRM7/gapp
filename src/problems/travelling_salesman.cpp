/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "travelling_salesman.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::problems
{
    auto TSP::invoke(const std::vector<PermutationGene>& chrom) const -> FitnessVector
    {
        double distance = 0.0;
        for (size_t i = 0; i < num_vars() - 1; i++)
        {
            distance += distance_matrix_[chrom[i]][chrom[i + 1]];
        }
        distance += distance_matrix_[chrom.front()][chrom.back()];

        return { -distance }; /* For maximization. */
    }

} // namespace genetic_algorithm::problems