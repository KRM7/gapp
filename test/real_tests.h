/* Benchmark/test functions for the real coded GA (RCGA). */

#ifndef REAL_TESTS_H
#define REAL_TESTS_H

#include <vector>
#include <utility>
#include <iostream>
#include <iomanip>

#include "../src/algorithms/real_ga.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/real.hpp"
#include "../src/mutation/real.hpp"
#include "../src/stop_condition/stop_condition.hpp"
#include "fitness_functions.h"
#include "utils.h"

using namespace std;
using namespace genetic_algorithm;

void realRastriginTest()
{
    Rastrigin rastriginFunction(10);

    pair<double, double> limit = { rastriginFunction.lbound(), rastriginFunction.ubound() };
    vector<pair<double, double>> limits(rastriginFunction.num_vars, limit);

    RCGA GA(rastriginFunction.num_vars, rastriginFunction, limits);

    GA.population_size(100);
    GA.selection_method(selection::single_objective::Roulette{});
    GA.crossover_method(crossover::real::SimulatedBinary{ limits, 0.6, 2.0 });
    GA.mutation_method(mutation::real::Gauss{ limits, 0.05 });
    GA.stop_condition(stopping::FitnessValue{ { -0.01 } });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 1000);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Rastrigin function is at (best is all " << rastriginFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << gene << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << rastriginFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void realRosenbrockTest()
{
    Rosenbrock rosenbrockFunction(10);

    pair<double, double> limit = { rosenbrockFunction.lbound(), rosenbrockFunction.ubound() };
    vector<pair<double, double>> limits(rosenbrockFunction.num_vars, limit);

    RCGA GA(rosenbrockFunction.num_vars, rosenbrockFunction, limits);

    GA.population_size(500);
    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::real::BLXa{ limits, 0.9 });
    GA.mutation_method(mutation::real::Uniform{ limits, 1.0 / rosenbrockFunction.num_vars });
    GA.stop_condition(stopping::FitnessEvals{ 500 * 1000 });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 2000);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Rosenbrock function is at (best is all " << rosenbrockFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << gene << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << rosenbrockFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void realSchwefelTest()
{
    Schwefel schwefelFunction(10);

    pair<double, double> limit = { schwefelFunction.lbound(), schwefelFunction.ubound() };
    vector<pair<double, double>> limits(schwefelFunction.num_vars, limit);

    RCGA GA(schwefelFunction.num_vars, schwefelFunction, limits);

    GA.population_size(500);
    GA.selection_method(selection::single_objective::Sigma{});
    GA.crossover_method(crossover::real::BLXa{ limits, 0.7 });
    GA.mutation_method(mutation::real::NonUniform{ limits, 1.0 / schwefelFunction.num_vars });
    GA.stop_condition(stopping::FitnessMeanStall{ 75, 0.01 });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 1000);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Schwefel function is at (best is all " << schwefelFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << gene << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << schwefelFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void realGriewankTest()
{
    Griewank griewankFunction(10);

    pair<double, double> limit = { griewankFunction.lbound(), griewankFunction.ubound() };
    vector<pair<double, double>> limits(griewankFunction.num_vars, limit);

    RCGA GA(griewankFunction.num_vars, griewankFunction, limits);

    GA.population_size(200);
    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::real::Wright{ limits, 0.85 });
    GA.mutation_method(mutation::real::Gauss{ limits, 0.05 });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 1500);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Griewank function is at (best is all " << griewankFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << gene << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << griewankFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void realAckleyTest()
{
    Ackley ackleyFunction(10);

    pair<double, double> limit = { ackleyFunction.lbound(), ackleyFunction.ubound() };
    vector<pair<double, double>> limits(ackleyFunction.num_vars, limit);

    RCGA GA(ackleyFunction.num_vars, ackleyFunction, limits);

    GA.population_size(200);
    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::real::Arithmetic{ limits, 0.85 });
    GA.mutation_method(mutation::real::Polynomial{ limits, 1.0 / ackleyFunction.num_vars, 60.0 });
    GA.stop_condition(stopping::FitnessBestStall{ 75, 0.002 });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 1000);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Ackley function is at (best is all " << ackleyFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        for (const auto& gene : sol.chromosome)
        {
            cout << gene << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << ackleyFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

#endif // !REAL_TESTS_H