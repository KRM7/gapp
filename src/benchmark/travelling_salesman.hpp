/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BENCHMARK_TSP_HPP
#define GA_BENCHMARK_TSP_HPP

#include "tsp_data/tsp_data.hpp"
#include "benchmark_function.hpp"
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <cstddef>

/* Traveling salesman problems to benchmark the PermutationGA. */
namespace genetic_algorithm::benchmark
{
    class TSP : public BenchmarkFunction<PermutationGene>
    {
    public:
        using Coords = std::array<double, 2>;
        using DistanceMatrix = std::vector<std::vector<double>>;

        template<size_t N>
        TSP(const std::array<Coords, N>& cities, double optimal_value) :
            BenchmarkFunction<PermutationGene>("TSP" + std::to_string(N), 1, N, Bounds{ 0, N - 1 }), distance_matrix_(N, std::vector(N, 0.0)), optimal_value_(optimal_value)
        {
            for (size_t i = 0; i < distance_matrix_.size(); i++)
                for (size_t j = 0; j < distance_matrix_.size(); j++)
                    distance_matrix_[i][j] = std::sqrt(std::pow(cities[i][0] - cities[j][0], 2) + std::pow(cities[i][1] - cities[j][1], 2));
        }

        double optimal_value() const noexcept { return optimal_value_; }

    private:
        std::vector<double> invoke(const std::vector<PermutationGene>& x) const override;

        DistanceMatrix distance_matrix_;
        double optimal_value_;
    };


    class TSP52 : public TSP
    {
    public:
        TSP52() : TSP(tsp52_coords, -7542.0) {}
    };


    class TSP76 : public TSP
    {
    public:
        TSP76() : TSP(tsp76_coords, -108159.0) {}
    };


    class TSP124 : public TSP
    {
    public:
        TSP124() : TSP(tsp124_coords, -59030.0) {}
    };


    class TSP152 : public TSP
    {
    public:
        TSP152() : TSP(tsp152_coords, -73682.0) {}
    };


    class TSP226 : public TSP
    {
    public:
        TSP226() : TSP(tsp226_coords, -80369.0) {}
    };


    class TSP299 : public TSP
    {
    public:
        TSP299() : TSP(tsp299_coords, -48191.0) {}
    };


    class TSP439 : public TSP
    {
    public:
        TSP439() : TSP(tsp439_coords, -107217.0) {}
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_TSP_HPP