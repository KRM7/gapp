/* Benchmark/test functions for the binary coded GA. */

#ifndef BINARY_TESTS_H
#define BINARY_TESTS_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include "../src/algorithms/binary_ga.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/binary.hpp"
#include "../src/mutation/binary.hpp"
#include "../src/stop_condition/stop_condition.hpp"
#include "fitness_functions.h"
#include "utils.h"

using namespace std;
using namespace genetic_algorithm;


void binaryRastriginTest()
{
    Rastrigin rastriginFunction(10);
    BinaryGA GA(400, rastriginFunction.num_vars * rastriginFunction.var_bits, rastriginFunction);
    
    GA.selection_method(selection::single_objective::Roulette{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.015 });
    GA.stop_condition(stopping::FitnessMeanStall{ 50, 0.005 });

    auto [sols, time_spent] = timed(&BinaryGA::run, GA, 1000);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Rastrigin function is at (best is all " << rastriginFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        vector<double> vars = convertToReals(sol.chromosome, rastriginFunction.var_bits, rastriginFunction.intval());
        for_each(vars.begin(), vars.end(), [&rastriginFunction](double& var) { var += rastriginFunction.lbound(); });

        for (const auto& var : vars)
        {
            cout << var << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << rastriginFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void binaryRosenbrockTest()
{
    Rosenbrock rosenbrockFunction(10);
    BinaryGA GA(300, rosenbrockFunction.num_vars * rosenbrockFunction.var_bits, rosenbrockFunction);

    GA.selection_method(selection::single_objective::Tournament{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.8 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });

    auto [sols, time_spent] = timed(&BinaryGA::run, GA, 1500);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Rosenbrock function is at (best is all " << rosenbrockFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        vector<double> vars = convertToReals(sol.chromosome, rosenbrockFunction.var_bits, rosenbrockFunction.intval());
        for_each(vars.begin(), vars.end(), [&rosenbrockFunction](double& var) { var += rosenbrockFunction.lbound(); });

        for (const auto& var : vars)
        {
            cout << var << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << rosenbrockFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void binarySchwefelTest()
{
    Schwefel schwefelFunction(10);
    BinaryGA GA(200, schwefelFunction.num_vars * schwefelFunction.var_bits, schwefelFunction);

    GA.selection_method(selection::single_objective::Rank{});
    GA.crossover_method(crossover::binary::Uniform{ 0.7 });
    GA.mutation_method(mutation::binary::Flip{ 0.01 });
    GA.stop_condition(stopping::FitnessEvals{ 200 * 1000 });

    auto [sols, time_spent] = timed(&BinaryGA::run, GA, 1500);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Schwefel function is at (best is all " << schwefelFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        vector<double> vars = convertToReals(sol.chromosome, schwefelFunction.var_bits, schwefelFunction.intval());
        for_each(vars.begin(), vars.end(), [&schwefelFunction](double& var) { var += schwefelFunction.lbound(); });

        for (const auto& var : vars)
        {
            cout << var << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << schwefelFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void binaryGriewankTest()
{
    Griewank griewankFunction(10);
    BinaryGA GA(250, griewankFunction.num_vars * griewankFunction.var_bits, griewankFunction);

    GA.selection_method(selection::single_objective::Sigma{});
    GA.crossover_method(crossover::binary::TwoPoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.04 });
    GA.stop_condition(stopping::FitnessValue{ { -0.1 } });

    auto [sols, time_spent] = timed(&BinaryGA::run, GA, 2500);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Griewank function is at (best is all " << griewankFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        vector<double> vars = convertToReals(sol.chromosome, griewankFunction.var_bits, griewankFunction.intval());
        for_each(vars.begin(), vars.end(), [&griewankFunction](double& var) { var += griewankFunction.lbound(); });

        for (const auto& var : vars)
        {
            cout << var << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << griewankFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

void binaryAckleyTest()
{
    Ackley ackleyFunction(10);
    BinaryGA GA(250, ackleyFunction.num_vars * ackleyFunction.var_bits, ackleyFunction);

    GA.selection_method(selection::single_objective::Boltzmann{});
    GA.crossover_method(crossover::binary::SinglePoint{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.04 });
    GA.stop_condition(stopping::FitnessBestStall{ 50, 0.002 });

    auto [sols, time_spent] = timed(&BinaryGA::run, GA, 2500);

    cout << setprecision(4);
    cout << "\n\nThe optimum of the Ackley function is at (best is all " << ackleyFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        vector<double> vars = convertToReals(sol.chromosome, ackleyFunction.var_bits, ackleyFunction.intval());
        for_each(vars.begin(), vars.end(), [&ackleyFunction](double& var) { var += ackleyFunction.lbound(); });

        for (const auto& var : vars)
        {
            cout << var << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << ackleyFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";
}

#endif // !BINARY_TESTS_H