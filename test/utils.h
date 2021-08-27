/* Functions used in the test functions. */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "../src/base_ga.h"

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

#endif // !UTILS_H