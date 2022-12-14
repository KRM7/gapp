/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "travelling_salesman.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::problems
{
    std::vector<double> TSP::invoke(const std::vector<PermutationGene>& x) const
    {
        double distance = 0.0;
        for (size_t i = 0; i < num_vars() - 1; i++)
        {
            distance += distance_matrix_[x[i]][x[i + 1]];
        }
        distance += distance_matrix_[x.front()][x.back()];

        return { -distance }; /* For maximization. */
    }

} // namespace genetic_algorithm::problems