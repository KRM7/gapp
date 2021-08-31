/* Benchmark/test functions for the binary coded GA. */

#ifndef BINARY_TESTS_H
#define BINARY_TESTS_H

#include <vector>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iomanip>

#include "../src/binary_ga.h"
#include "fitness_functions.h"
#include "utils.h"

using namespace std;
using namespace genetic_algorithm;


void binaryRastriginTest()
{
    /* Init GA. */
    Rastrigin rastriginFunction(10);
    BinaryGA GA(rastriginFunction.num_vars * rastriginFunction.var_bits, rastriginFunction);
    
    /* Set some optional parameters. */
    GA.population_size(400);
    GA.crossover_rate(0.75);
    GA.mutation_rate(0.015);
    GA.selection_method(BinaryGA::SogaSelection::roulette);
    GA.crossover_method(BinaryGA::CrossoverMethod::two_point);

    GA.max_gen(1500);
    GA.stop_condition(BinaryGA::StopCondition::fitness_mean_stall);
    GA.stall_gen_count(50);
    GA.stall_threshold(0.005);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << setprecision(4);
    cout << "\n\nThe optimum of the Rastrigin function is at (best is all " << rastriginFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        /* Decode sol. */
        vector<double> vars = convertToReals(sol.chromosome, rastriginFunction.var_bits, rastriginFunction.intval());
        for_each(vars.begin(), vars.end(), [&rastriginFunction](double& var) { var += rastriginFunction.lbound(); });
        /* Print. */
        for (const auto& var : vars)
        {
            cout << var << "  ";
        }
        cout << "\n";
    }
    cout << "Fitness value: " << sols[0].fitness[0] << " (best is " << rastriginFunction.optimal_value() << ")\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << time_spent << " s\n\n";

    /* Print the GA stats. */
    //displayStats(GA.soga_history());
}

void binaryRosenbrockTest()
{
    /* Init GA. */
    Rosenbrock rosenbrockFunction(10);
    BinaryGA GA(rosenbrockFunction.num_vars * rosenbrockFunction.var_bits, rosenbrockFunction);

    /* Set some optional parameters. */
    GA.population_size(300);
    GA.crossover_rate(0.8);
    GA.mutation_rate(0.01);
    GA.selection_method(BinaryGA::SogaSelection::tournament);
    GA.crossover_method(BinaryGA::CrossoverMethod::n_point);
    GA.num_crossover_points(2);
    GA.max_gen(1500);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << setprecision(4);
    cout << "\n\nThe optimum of the Rosenbrock function is at (best is all " << rosenbrockFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        /* Decode sol */
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

    /* Print the GA stats. */
    //displayStats(GA.soga_history());
}

void binarySchwefelTest()
{
    /* Init GA. */
    Schwefel schwefelFunction(10);
    BinaryGA GA(schwefelFunction.num_vars * schwefelFunction.var_bits, schwefelFunction);

    /* Set some optional parameters. */
    GA.population_size(200);
    GA.crossover_rate(0.7);
    GA.mutation_rate(0.01);
    GA.selection_method(BinaryGA::SogaSelection::rank);
    GA.crossover_method(BinaryGA::CrossoverMethod::uniform);
    
    GA.max_gen(1500);
    GA.stop_condition(BinaryGA::StopCondition::fitness_evals);
    GA.max_fitness_evals(200 * 650);

    /* Run GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << setprecision(4);
    cout << "\n\nThe optimum of the Schwefel function is at (best is all " << schwefelFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        /* Decode sol. */
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

    /* Print GA stats. */
    //displayStats(GA.soga_history());
}

void binaryGriewankTest()
{
    /* Init GA. */
    Griewank griewankFunction(10);
    BinaryGA GA(griewankFunction.num_vars * griewankFunction.var_bits, griewankFunction);

    /* Set some optional parameters. */
    GA.population_size(250);
    GA.crossover_rate(0.75);
    GA.mutation_rate(0.04);
    GA.selection_method(BinaryGA::SogaSelection::sigma);
    GA.crossover_method(BinaryGA::CrossoverMethod::n_point);
    GA.num_crossover_points(3);

    GA.max_gen(2500);
    GA.stop_condition(BinaryGA::StopCondition::fitness_value);
    GA.fitness_threshold({ -0.1 });

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << setprecision(4);
    cout << "\n\nThe optimum of the Griewank function is at (best is all " << griewankFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        /* Decode sol. */
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

    /* Print GA stats. */
    //displayStats(GA.soga_history());
}

void binaryAckleyTest()
{
    /* Init GA. */
    Ackley ackleyFunction(10);
    BinaryGA GA(ackleyFunction.num_vars * ackleyFunction.var_bits, ackleyFunction);

    /* Set some optional parameters. */
    GA.population_size(250);
    GA.crossover_rate(0.75);
    GA.mutation_rate(0.04);
    GA.selection_method(BinaryGA::SogaSelection::boltzmann);
    GA.crossover_method(BinaryGA::CrossoverMethod::single_point);
    
    GA.max_gen(2500);
    GA.stop_condition(BinaryGA::StopCondition::fitness_best_stall);
    GA.stall_gen_count(50);
    GA.stall_threshold(0.002);

    /* Run the GA with a timer. */
    auto tbegin = chrono::high_resolution_clock::now();
    auto sols = GA.run();
    auto tend = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(tend - tbegin).count();
    double time_spent = double(duration) / 1E+6;

    /* Print the results. */
    cout << setprecision(4);
    cout << "\n\nThe optimum of the Ackley function is at (best is all " << ackleyFunction.optimal_x() << "): \n";
    for (const auto& sol : sols)
    {
        /* Decode sol. */
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

    /* Print GA stats. */
    //displayStats(GA.soga_history());
}

#endif // !BINARY_TESTS_H