#ifndef TIME_TEST_H
#define TIME_TEST_H

#include "../src/algorithms/binary_ga.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/binary.hpp"
#include "../src/mutation/binary.hpp"
#include "utils.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace genetic_algorithm;

void timeGA(size_t num_runs = 50)
{
    BinaryGA GA(100, 100,
    [](const vector<char>& chrom)
    {
        return vector(1, double(count(chrom.begin(), chrom.end(), 1)));
    });

    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.6 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });

    cout << setprecision(4);
    double running_mean_time = 0.0;
    for (size_t i = 1; i <= num_runs; i++)
    {
        auto [_, time_spent] = timed(&BinaryGA::run, GA, 1000);

        running_mean_time += (time_spent - running_mean_time) / i;
        cout << "Time taken: " << running_mean_time << " s \r";
    }
}

#endif // !TIME_TEST_H