/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/**
* Implementations of travelling salesman problems that can be used to test the
* single-objective permutational algorithms. \n
* All of the problems are implemented for maximization, so they return negative
* distances.
*/

#ifndef GA_PROBLEMS_TSP_HPP
#define GA_PROBLEMS_TSP_HPP

#include "tsp_data/tsp_data.hpp"
#include "benchmark_function.hpp"
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <cstddef>

namespace genetic_algorithm::problems
{
    /** ... */
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

    /**
    * Travelling salesman problem with 52 nodes (Berlin52) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP52 : public TSP
    {
    public:
        /** Default constructor. */
        TSP52() : TSP(tsp52_coords, -7542.0) {}
    };


    /**
    * Travelling salesman problem with 76 nodes (Padberg/Rinaldi's 76 city problem) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP76 : public TSP
    {
    public:
        /** Default constructor. */
        TSP76() : TSP(tsp76_coords, -108159.0) {}
    };

    /**
    * Travelling salesman problem with 124 nodes (Padberg/Rinaldi's 124 city problem) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP124 : public TSP
    {
    public:
        /** Default constructor. */
        TSP124() : TSP(tsp124_coords, -59030.0) {}
    };
    
    /**
    * Travelling salesman problem with 152 nodes (Padberg/Rinaldi's 152 city problem) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP152 : public TSP
    {
    public:
        /** Default constructor. */
        TSP152() : TSP(tsp152_coords, -73682.0) {}
    };

    /**
    * Travelling salesman problem with 226 nodes (Padberg/Rinaldi's 226 city problem) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP226 : public TSP
    {
    public:
        /** Default constructor. */
        TSP226() : TSP(tsp226_coords, -80369.0) {}
    };

    /**
    * Travelling salesman problem with 299 nodes (Padberg/Rinaldi's 299 city problem) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP299 : public TSP
    {
    public:
        /** Default constructor. */
        TSP299() : TSP(tsp299_coords, -48191.0) {}
    };

    /**
    * Travelling salesman problem with 439 nodes (Padberg/Rinaldi's 439 city problem) for testing the PermutationGA. \n
    * Implemented for maximization (returns negative distances).
    */
    class TSP439 : public TSP
    {
    public:
        /** Default constructor. */
        TSP439() : TSP(tsp439_coords, -107217.0) {}
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_TSP_HPP