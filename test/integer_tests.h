#ifndef INTEGER_TEST_H
#define INTEGER_TEST_H

#include <vector>
#include <chrono>
#include <iostream>

#include "../src/algorithms/integer_ga.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/integer.hpp"
#include "../src/mutation/integer.hpp"
#include "fitness_functions.h"
#include "utils.h"

using namespace std;
using namespace genetic_algorithm;

void integerTest1()
{
    MatchString match("HELLO WORLD!");

    IntegerGA GA(100, match.num_vars(), match, 96);

    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::integer::TwoPoint{ 0.85 });
    GA.mutation_method(mutation::integer::Uniform{ 0.01 });

    auto [sols, time_spent] = timed(&IntegerGA::run, GA, 500);

    cout << "\n\nThe best strings found are (expected: " << match.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << char(gene + 32);
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << match.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void integerTest2()
{
    MatchString match("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt.");

    IntegerGA GA(250, match.num_vars(), match, 96);

    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::integer::Uniform{ 0.8 });
    GA.mutation_method(mutation::integer::Uniform{ 5.0 / 250 });

    auto [sols, time_spent] = timed(&IntegerGA::run, GA, 1000);

    cout << "\n\nThe best strings found are (expected: " << match.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << char(gene + 32);
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << match.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

#endif // !INTEGER_TEST_H