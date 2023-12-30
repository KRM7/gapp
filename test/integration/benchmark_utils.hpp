/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_TEST_BENCHMARK_UTILS_HPP
#define GA_TEST_BENCHMARK_UTILS_HPP

#include "gapp.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <utility>
#include <atomic>
#include <type_traits>

using namespace gapp;
using namespace gapp::problems;

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

    return std::make_pair(std::move(result), time_spent);
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

template<typename T, typename F>
void benchmarkSoga(GA<T>& ga, size_t max_gen, F fitness_func)
{
    using RunFn = std::conditional_t<!is_bounded<T>,
        Candidates<T>(GA<T>::*)(F, size_t, Population<T>),
        Candidates<T>(GA<T>::*)(F, BoundsVector<T>, size_t, Population<T>)>;

    auto [sols, time_spent] = [&]
    {
        if constexpr (!is_bounded<T>)
            return invoke_timed(static_cast<RunFn>(&GA<T>::template solve<F>), ga, fitness_func, max_gen, Population<T>{});
        else
            return invoke_timed(static_cast<RunFn>(&GA<T>::template solve<F>), ga, fitness_func, fitness_func.bounds(), max_gen, Population<T>{});
    }();

    std::string algo = std::is_floating_point_v<T> ? "RCGA" : "BinaryGA";

    std::cout << "Function: " << fitness_func.name() << ", " << algo
              << "\nNumber of optimal sols: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << fitness_func.optimal_value()[0] << ")"
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

template<typename T, typename F>
void benchmarkMoga(GA<T>& ga, size_t max_gen, const std::string& ga_name, F fitness_func)
{
    using RunFn = std::conditional_t<!is_bounded<T>,
        Candidates<T>(GA<T>::*)(F, size_t, Population<T>),
        Candidates<T>(GA<T>::*)(F, BoundsVector<T>, size_t, Population<T>)>;

    auto [sols, time_spent] = [&]
    {
        if constexpr (!is_bounded<T>)
            return invoke_timed(static_cast<RunFn>(&GA<T>::template solve<F>), ga, fitness_func, max_gen, Population<T>{});
        else
            return invoke_timed(static_cast<RunFn>(&GA<T>::template solve<F>), ga, fitness_func, fitness_func.bounds(), max_gen, Population<T>{});
    }();

    std::string problem = fitness_func.name();

    std::cout << "Function: " << problem
              << "\nAlgorithm: " << ga_name
              << "\nNumber of optimal sols: " << sols.size()
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << time_spent << "s\n\n";

    std::string name = problem.substr(0, problem.find(','));
    std::string pop_file = "../tools/mo_results/" + ga_name + "_" + name + "_last.txt";
    std::string sol_file = "../tools/mo_results/" + ga_name + "_" + name + "_sols.txt";

    std::ofstream flast(pop_file);
    std::ofstream fsols(sol_file);

    writePopulationToFile(ga.population(), flast);
    writePopulationToFile(sols, fsols);
}

template<typename F>
void benchmarkTSP(PermutationGA& ga, size_t max_gen, F fitness_func)
{
    using RunFn = Candidates<PermutationGene>(PermutationGA::*)(F, size_t, Population<PermutationGene>);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&PermutationGA::template solve<F>), ga, fitness_func, max_gen, Population<PermutationGene>{});

    std::cout << "Function: " << fitness_func.name()
              << "\nNumber of optimal sols: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << fitness_func.optimal_value()[0] << ")"
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

template<typename F>
void benchmarkInt(IntegerGA& ga, size_t max_gen, F fitness_func)
{
    using RunFn = Candidates<IntegerGene>(IntegerGA::*)(F, BoundsVector<IntegerGene>, size_t, Population<IntegerGene>);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&IntegerGA::template solve<F>), ga, fitness_func, fitness_func.bounds(), max_gen, Population<IntegerGene>{});

    std::cout << "Function: " << fitness_func.name()
              << "\nNumber of optimal sols: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << fitness_func.optimal_value()[0] << ")"
              << "\nNumber of objective function evals: " << ga.num_fitness_evals()
              << " (instead of: " << max_gen * ga.population_size() << ")"
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

#endif // !GA_TEST_BENCHMARK_UTILS_HPP