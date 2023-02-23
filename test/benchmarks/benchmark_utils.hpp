/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_UTILS_HPP
#define GA_TEST_BENCHMARK_UTILS_HPP

#include "genetic_algorithm.hpp"
#include "problems/benchmark_function.hpp"
#include "problems/single_objective.hpp"
#include "problems/multi_objective.hpp"
#include "problems/travelling_salesman.hpp"
#include "problems/integer.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <utility>
#include <atomic>
#include <type_traits>

using namespace genetic_algorithm;
using namespace genetic_algorithm::problems;

template<typename F, typename... Args>
auto invoke_timed(F&& f, Args&&... args)
{
    const auto tbegin = std::chrono::high_resolution_clock::now();
    std::atomic_signal_fence(std::memory_order::seq_cst);
    auto result = std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    std::atomic_signal_fence(std::memory_order::seq_cst);
    const auto tend = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tbegin).count();
    double time_spent = duration / 1000.0;

    return std::make_pair(result, time_spent);
}

template<typename T>
void writePopulationToFile(const std::vector<T>& sols, std::ostream& os)
{
    for (const auto& sol : sols)
    {
        for (const auto& f : sol.fitness) { os << f << "\t"; }
        os << "\n";
    }
}

template<typename T>
void benchmarkSoga(GA<T>& ga, size_t max_gen, const BenchmarkFunction<RealGene>& fitness_func)
{
    using RunFn = Candidates<T>(GA<T>::*)(size_t);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&GA<T>::run), ga, max_gen);

    std::string algo = std::is_floating_point_v<T> ? "RCGA" : "BinaryGA";

    std::cout << "Function: " << fitness_func.name() << ", " << algo
              << "\nNumber of optimal sols: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << fitness_func.optimal_value()[0] << ")"
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

template<typename T>
void benchmarkMoga(GA<T>& ga, size_t max_gen, const std::string& ga_name, const BenchmarkFunction<T>& fitness_func)
{
    using RunFn = Candidates<T>(GA<T>::*)(size_t);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&GA<T>::run), ga, max_gen);

    std::string problem = fitness_func.name();

    std::cout << "Function: " << problem
              << "\nAlgorithm: " << ga_name
              << "\nNumber of optimal sols: " << sols.size()
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << time_spent << "s\n\n";

    std::string name = problem.substr(0, problem.find(','));
    std::string pop_file = "mo_results/" + ga_name + "_" + name + "_last.txt";
    std::string sol_file = "mo_results/" + ga_name + "_" + name + "_sols.txt";

    std::ofstream flast(pop_file);
    std::ofstream fsols(sol_file);

    writePopulationToFile(ga.population(), flast);
    writePopulationToFile(sols, fsols);
}

inline void benchmarkTSP(PermutationGA& ga, size_t max_gen, const TSP& fitness_func)
{
    using RunFn = Population<PermutationGene>(PermutationGA::*)(size_t);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&PermutationGA::run), ga, max_gen);

    std::cout << "Function: " << fitness_func.name()
              << "\nNumber of optimal sols: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << fitness_func.optimal_value()[0] << ")"
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

inline void benchmarkInt(IntegerGA& ga, size_t max_gen, const StringFinder& fitness_func)
{
    using RunFn = Candidates<IntegerGene>(IntegerGA::*)(size_t);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&IntegerGA::run), ga, max_gen);

    std::cout << "Function: " << fitness_func.name()
              << "\nNumber of optimal sols: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << fitness_func.optimal_value()[0] << ")"
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

#endif // !GA_TEST_BENCHMARK_UTILS_HPP