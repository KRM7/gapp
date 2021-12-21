/* Functions for measuring the speed of the GA. */

#ifndef TIME_TEST_H
#define TIME_TEST_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <chrono>

#include "../src/binary_ga.h"
#include "../src/crossover/binary.hpp"
#include "../src/mutation/binary.hpp"

using namespace std;
using namespace genetic_algorithm;

void timeGA(size_t num_runs = 50)
{
    /* Init GA. */
    auto fitness_func = [](const vector<char>& chrom) -> vector<double>
    { 
        return vector<double>(1, double(count(chrom.begin(), chrom.end(), 1))); 
    };

    BinaryGA GA(100, fitness_func);

    /* Set some optional parameters. */
    GA.population_size(100);
    GA.selection_method(BinaryGA::SogaSelection::tournament);
    GA.crossover_method(crossover::binary::TwoPoint{ 0.6 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.max_gen(1000);

    /* Run the GA with a timer. */
    cout << setprecision(4);
    double running_mean_time = 0.0;
    for (size_t i = 1; i <= num_runs; i++)
    {
        auto tbegin = chrono::high_resolution_clock::now();
        GA.run();
        auto tend = chrono::high_resolution_clock::now();

        auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
        double time_spent = double(duration) / 1E+6;

        running_mean_time += (time_spent - running_mean_time) / i;
        cout << "Time taken: " << running_mean_time << " s \r";
    }
}

#endif // !TIME_TEST_H