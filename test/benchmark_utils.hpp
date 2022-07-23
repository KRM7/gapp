/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef BENCHMARK_UTILS_HPP
#define BENCHMARK_UTILS_HPP

#include "../src/genetic_algorithm.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <format>
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
    auto tbegin = std::chrono::high_resolution_clock::now();
    std::atomic_signal_fence(std::memory_order::seq_cst);
    auto result = std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    std::atomic_signal_fence(std::memory_order::seq_cst);
    auto tend = std::chrono::high_resolution_clock::now();

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
    size_t var_count = binary_chrom.size() / bits_per_var;

    std::vector<double> vars(var_count);

    for (size_t i = 0; i < vars.size(); i++)
    {
        auto first = binary_chrom.begin() + i * bits_per_var;
        auto last = binary_chrom.begin() + (i + 1) * bits_per_var;

        double val = std::accumulate(first, last, 0.0,
        [](double acc, char bit) noexcept
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
    auto [sols, time_spent] = invoke_timed(&GA::run, ga, max_gen);

    std::cout 
    << std::format("\n\nOptimum found for the {} is (actual best is {}):\n", problem_name, fitness_func.optimal_x());

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

    std::cout
    << std::format("The number of optimal solutions found: {}\n"
                   "Best fitness found: {:.4g} (best possible is {:.4g})\n"
                   "Number of objective function evals performed: {}\n"
                   "Time taken: {:.4f} s\n\n",
                   sols.size(),
                   sols[0].fitness[0], fitness_func.optimal_value(),
                   ga.num_fitness_evals(),
                   time_spent);
}

template<genetic_algorithm::GeneticAlgorithmType GA>
void benchmarkMoga(GA& ga, size_t max_gen, const std::string& ga_name, const std::string& problem_name)
{
    auto [sols, time_spent] = invoke_timed(&GA::run, ga, max_gen);

    std::string msg = "\n\nOptimal solutions found for the {} problem with the {}: {}\n"
                      "Number of fitness function evaluations: {}\n"
                      "Time taken: {:.4f} s\n\n";

    std::cout << std::format(msg, problem_name, ga_name, sols.size(), ga.num_fitness_evals(), time_spent);

    std::string sols_path = std::format("test/mo_results/{}_{}_{}.txt", ga_name, problem_name, "{}");

    std::ofstream flast(std::format(sols_path, "last"));
    std::ofstream fsols(std::format(sols_path, "sols"));

    writePopulationToFile(ga.population(), flast);
    writePopulationToFile(sols, fsols);
}

#endif // !BENCHMARK_UTILS_HPP