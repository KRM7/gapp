/* Benchmarks for the NSGA-II algorithm. */

#ifndef NSGA2_BENCHMARK_HPP
#define NSGA2_BENCHMARK_HPP

#include "../src/encoding/real.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/real.hpp"
#include "../src/mutation/real.hpp"
#include "fitness_functions.h"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;

void nsga2_kur()
{
    KUR kur(3);

    RCGA GA(kur.num_vars, kur, kur.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / kur.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "KUR");
}

void nsga2_zdt2()
{
    ZDT2 zdt2(30);
    
    RCGA GA(zdt2.num_vars, zdt2, zdt2.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / zdt2.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "ZDT2");
}

void nsga2_zdt3()
{
    ZDT3 zdt3(30);

    RCGA GA(zdt3.num_vars, zdt3, zdt3.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / zdt3.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "ZDT3");
}

void nsga2_zdt6()
{
    ZDT6 zdt6(10);
    
    RCGA GA(zdt6.num_vars, zdt6, zdt6.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / zdt6.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA2", "ZDT6");
}

void nsga2_dtlz1()
{
    DTLZ1 dtlz1(7, 3);

    RCGA GA(dtlz1.num_vars, dtlz1, dtlz1.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / dtlz1.bounds().size() });

    benchmarkMoga(GA, 1500, "NSGA2", "DTLZ1");
}

void nsga2_dtlz2()
{
    DTLZ2 dtlz2(12, 3);

    RCGA GA(dtlz2.num_vars, dtlz2, dtlz2.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA2{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / dtlz2.bounds().size() });

    benchmarkMoga(GA, 1500, "NSGA2", "DTLZ2");
}


#endif // !NSGA2_BENCHMARK_HPP