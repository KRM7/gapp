/* Tests for the NSGA-II algorithm. */

#ifndef NSGA2_TESTS_H
#define NSGA2_TESTS_H

#include <vector>
#include <utility>
#include <chrono>
#include <numeric>
#include <cmath>
#include <cassert>
#include <iostream>

#include "../src/algorithms/real_ga.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/real.hpp"
#include "../src/mutation/real.hpp"
#include "fitness_functions.h"
#include "utils.h"

using namespace std;
using namespace genetic_algorithm;


void nsga2KurTest()
{
    KUR kurFunction(3);
    pair<double, double> limit = { kurFunction.lbound(), kurFunction.ubound() };
    vector<pair<double, double>> limits(kurFunction.num_vars, limit);

    RCGA GA(kurFunction.num_vars, kurFunction, limits);

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ limits, 1.0 / limits.size() });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 250);

    cout << "\n\nNumber of optimal solutions found for the KUR problem with the NSGA-II: " << sols.size() << "\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << setprecision(4) << time_spent << " s\n\n";

    //writeResultsToFile(GA.population(), "test/mo_results/nsga2_kur_last.txt");
    writeResultsToFile(sols, "test/mo_results/nsga2_kur_sols.txt");
}

void nsga2Zdt2Test()
{
    ZDT2 zdt2Function(30);
    pair<double, double> limit = { zdt2Function.lbound(), zdt2Function.ubound() };
    vector<pair<double, double>> limits(zdt2Function.num_vars, limit);
    
    RCGA GA(zdt2Function.num_vars, zdt2Function, limits);

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ limits, 1.0 / limits.size() });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 250);

    cout << "\n\nNumber of optimal solutions found for the ZDT2 problem with the NSGA-II: " << sols.size() << "\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << setprecision(4) << time_spent << " s\n\n";

    //writeResultsToFile(GA.population(), "test/mo_results/nsga2_zdt2_last.txt");
    writeResultsToFile(sols, "test/mo_results/nsga2_zdt2_sols.txt");
}

void nsga2Zdt3Test()
{
    ZDT3 zdt3Function(30);
    pair<double, double> limit = { zdt3Function.lbound(), zdt3Function.ubound() };
    vector<pair<double, double>> limits(zdt3Function.num_vars, limit);

    RCGA GA(zdt3Function.num_vars, zdt3Function, limits);

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ limits, 1.0 / limits.size() });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 250);

    cout << "\n\nNumber of optimal solutions found for the ZDT3 problem with the NSGA-II: " << sols.size() << "\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << setprecision(4) << time_spent << " s\n\n";

    //writeResultsToFile(GA.population(), "test/mo_results/nsga2_zdt3_last.txt");
    writeResultsToFile(sols, "test/mo_results/nsga2_zdt3_sols.txt");
}

void nsga2Zdt6Test()
{
    ZDT6 zdt6Function(10);
    pair<double, double> limit = { zdt6Function.lbound(), zdt6Function.ubound() };
    vector<pair<double, double>> limits(zdt6Function.num_vars, limit);
    
    RCGA GA(zdt6Function.num_vars, zdt6Function, limits);

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ limits, 1.0 / limits.size() });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 250);

    cout << "\n\nNumber of optimal solutions found for the ZDT6 problem with the NSGA-II: " << sols.size() << "\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << setprecision(4) << time_spent << " s\n\n";

    //writeResultsToFile(GA.population(), "test/mo_results/nsga2_zdt6_last.txt");
    writeResultsToFile(sols, "test/mo_results/nsga2_zdt6_sols.txt");
}

void nsga2Dtlz1Test()
{
    DTLZ1 dtlz1Function(7, 3);
    pair<double, double> limit = { dtlz1Function.lbound(), dtlz1Function.ubound() };
    vector<pair<double, double>> limits(dtlz1Function.num_vars, limit);

    RCGA GA(dtlz1Function.num_vars, dtlz1Function, limits);

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ limits, 1.0 / limits.size() });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 1500);

    cout << "\n\nNumber of optimal solutions found for the DTLZ1 problem with the NSGA-II: " << sols.size() << "\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << setprecision(4) << time_spent << " s\n\n";

    //writeResultsToFile(GA.population(), "test/mo_results/nsga2_dtlz1_last.txt");
    writeResultsToFile(sols, "test/mo_results/nsga2_dtlz1_sols.txt");
}

void nsga2Dtlz2Test()
{
    DTLZ2 dtlz2Function(12, 3);
    pair<double, double> limit = { dtlz2Function.lbound(), dtlz2Function.ubound() };
    vector<pair<double, double>> limits(dtlz2Function.num_vars, limit);

    RCGA GA(dtlz2Function.num_vars, dtlz2Function, limits);

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ limits, 1.0 / limits.size() });

    auto [sols, time_spent] = timed(&RCGA::run, GA, 1500);

    cout << "\n\nNumber of optimal solutions found for the DTLZ2 problem with the NSGA-II: " << sols.size() << "\n";
    cout << "Number of fitness evals: " << GA.num_fitness_evals() << "\n";
    cout << "Time taken: " << setprecision(4) << time_spent << " s\n\n";

    //writeResultsToFile(GA.population(), "test/mo_results/nsga2_dtlz2_last.txt");
    writeResultsToFile(sols, "test/mo_results/nsga2_dtlz2_sols.txt");
}

#endif // !NSGA2_TESTS_H