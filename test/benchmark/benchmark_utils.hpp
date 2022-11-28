/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef BENCHMARK_UTILS_HPP
#define BENCHMARK_UTILS_HPP

#include "../src/genetic_algorithm.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <utility>
#include <functional>
#include <algorithm>
#include <atomic>
#include <type_traits>


template<typename F, typename... Args>
requires std::invocable<F, Args...>
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
        for (const auto& f : sol.fitness)
        {
            os << f << "\t";
        }
        os << "\n";
    }
}

auto convertToReals(const std::vector<char>& binary_chrom, size_t bits_per_var, double interval_len, double lower_bound)
{
    const size_t var_count = binary_chrom.size() / bits_per_var;

    std::vector<double> vars(var_count);

    for (size_t i = 0; i < vars.size(); i++)
    {
        const auto first = binary_chrom.begin() + i * bits_per_var;
        const auto last = binary_chrom.begin() + (i + 1) * bits_per_var;

        const double val = std::accumulate(first, last, 0.0, [](double acc, char bit) noexcept
        {
            return (acc * 2) + bit;
        });

        vars[i] = val / (std::pow(2.0, bits_per_var) - 1);
        vars[i] *= interval_len;
        vars[i] += lower_bound;
    }

    return vars;
}

template<typename T>
void printSol(const std::vector<T>& chrom)
{
    std::cout << std::fixed << std::setw(6);
    for (const auto& gene : chrom)
    {
        std::cout << gene << "  ";
    }
    std::cout << "\n";
}

template<genetic_algorithm::GeneticAlgorithmType GA, typename F>
void benchmarkSoga(GA& ga, size_t max_gen, const F& fitness_func, const std::string& problem_name)
{
    using RunFn = genetic_algorithm::Candidates<typename GA::GeneType>(GA::*)(size_t);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&GA::run), ga, max_gen);

    std::cout << "\n\nOptimum found for the " << problem_name << " is (actual best is " << fitness_func.optimal_x() << "):\n";

    for (const auto& sol : sols)
    {
        if constexpr (std::is_same_v<GA, genetic_algorithm::BinaryGA>)
        {
            auto real_chrom = convertToReals(sol.chromosome, fitness_func.var_bits, fitness_func.intval(), fitness_func.lbound());
            printSol(real_chrom);
        }
        else if constexpr (!std::is_same_v<GA, genetic_algorithm::PermutationGA>)
        {
            printSol(sol.chromosome);
        }
    }

    std::cout << "The number of optimal solutions found: " << sols.size()
              << "\nBest fitness found: " << std::defaultfloat << std::setprecision(4) << sols[0].fitness[0]
              << " (best possible is " << std::defaultfloat << fitness_func.optimal_value()
              << "\nNumber of objective function evals performed: " << ga.num_fitness_evals()
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << "s\n\n";
}

template<genetic_algorithm::GeneticAlgorithmType GA>
void benchmarkMoga(GA& ga, size_t max_gen, const std::string& ga_name, const std::string& problem_name)
{
    using RunFn = genetic_algorithm::Candidates<typename GA::GeneType>(GA::*)(size_t);

    auto [sols, time_spent] = invoke_timed(static_cast<RunFn>(&GA::run), ga, max_gen);

    std::cout << "\n\nOptimal solutions found for the " << problem_name << " problem with the " << ga_name << ": " << sols.size() << "\n"
              << "\nNumber of fitness function evaluations: " << ga.num_fitness_evals()
              << "\nTime taken: " << std::fixed << std::setprecision(4) << time_spent << " s\n\n";

    std::string pop_path = "test/mo_results/" + ga_name + "_" + problem_name + "_last.txt";
    std::string sol_path = "test/mo_results/" + ga_name + "_" + problem_name + "_sols.txt";

    std::ofstream flast(pop_path);
    std::ofstream fsols(sol_path);

    writePopulationToFile(ga.population(), flast);
    writePopulationToFile(sols, fsols);
}

#endif // !BENCHMARK_UTILS_HPP