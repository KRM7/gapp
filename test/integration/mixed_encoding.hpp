/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_TEST_MIXED_ENCODING_HPP
#define GAPP_TEST_MIXED_ENCODING_HPP

#include "gapp.hpp"
#include <span>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace gapp;

class VRP : public FitnessFunctionBase<MixedGene<PermutationGene, IntegerGene>>
{
public:
    using Coords = std::array<double, 2>;

    VRP(std::span<const Coords> cities, size_t agents = 2) :
        FitnessFunctionBase({ cities.size(), cities.size() }),
        cities_(cities.begin(), cities.end()),
        agents_(agents)
    {}

private:
    FitnessVector invoke(const Candidate<MixedGene<PermutationGene, IntegerGene>>& sol) const override
    {
        small_vector<double> distances(agents_, 0.0);
        small_vector<Coords> positions(agents_, { 0.0, 0.0 });

        for (size_t i = 0; i < sol.chrom_len<PermutationGene>(); i++)
        {
            const size_t city = static_cast<size_t>(sol.chrom<PermutationGene>()[i]);
            const size_t agent = static_cast<size_t>(sol.chrom<IntegerGene>()[city]);

            const Coords old_pos = positions[agent];
            const Coords new_pos = cities_[city];

            distances[agent] += std::hypot(new_pos[0] - old_pos[0], new_pos[1] - old_pos[1]);
            positions[agent] = new_pos;
        }

        const double max_distance = *std::max_element(distances.begin(), distances.end());

        return { -max_distance };
    }

    std::vector<Coords> cities_;
    size_t agents_;
};


template<auto Coords, size_t Agents = 2>
inline void mixed_encoding()
{
    MixedGA<PermutationGene, IntegerGene> ga(500);

    ga.crossover_method(crossover::Mixed{ crossover::perm::Order2{ 0.9 }, crossover::integer::Uniform{} });
    ga.mutation_method(mutation::Mixed{ mutation::perm::Inversion{ 0.5 }, mutation::integer::Uniform{ 0.2 } });

    const auto solutions = ga.solve(VRP{ Coords, Agents }, Bounds<IntegerGene>{ 0, Agents - 1 }, 1500);

    std::cout << "\nMixed encoded problem with size=" << Coords.size() << ", agents=" << Agents << ":"
        << "\nNumber of optimal sols: " << solutions.size()
        << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << solutions[0].fitness[0] << "\n";
}

#endif // !GAPP_TEST_MIXED_ENCODING_HPP
