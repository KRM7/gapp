#ifndef UTILS_H
#define UTILS_H

#include "../src/algorithms/ga_base.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <utility>
#include <functional>
#include <atomic>

using namespace genetic_algorithm;

/* Prints the statistics of a genetic algorithm. */
void displayStats(auto stats)
{
    std::cout << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < stats.fitness_mean.size(); i++)
    {
        if (i % 20 == 0)
        {
            std::cout << "****************************************" << "\n"
                      << " gen |   avgf    |    maxf   |   fSD " << "\n"
                      << "****************************************" << "\n";
        }
        std::cout << i + 1 << " | " << stats.fitness_mean[i] << " | " << stats.fitness_max[i] << " | " << stats.fitness_sd[i] << "\n";
    }
}

/* Write the fitness values of the solutions to a file for plotting. */
void writeResultsToFile(auto sols, std::string fname)
{
    std::ofstream file;

    file.open(fname);
    for (const auto& sol : sols)
    {
        for (const auto& f : sol.fitness)
        {
            file << f << "\t";
        }
        file << "\n";
    }
    file.close();
}

template<typename F, typename... Args>
auto timed(F&& f, Args&&... args)
{
    auto tbegin = std::chrono::high_resolution_clock::now();
    std::atomic_signal_fence(std::memory_order_seq_cst);
    auto result = std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    std::atomic_signal_fence(std::memory_order_seq_cst);
    auto tend = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    return std::make_pair(result, time_spent);
}

#endif // !UTILS_H