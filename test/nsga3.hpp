/* Tests for the NSGA-III algorithm. */

#ifndef NSGA3_BENCHMARK_HPP
#define NSGA3_BENCHMARK_HPP

#include "../src/encoding/real.hpp"
#include "../src/selection/selection.hpp"
#include "../src/crossover/real.hpp"
#include "../src/mutation/real.hpp"
#include "fitness_functions.h"
#include "benchmark_utils.hpp"

using namespace genetic_algorithm;

void nsga3_kur()
{
    KUR kur(3);

    RCGA GA(kur.num_vars, kur, kur.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / kur.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA3", "KUR");
}

void nsga3_zdt2()
{
    ZDT2 zdt2(30);

    RCGA GA(zdt2.num_vars, zdt2, zdt2.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / zdt2.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA3", "ZDT2");
}

void nsga3_zdt3()
{
    ZDT3 zdt3(30);

    RCGA GA(zdt3.num_vars, zdt3, zdt3.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / zdt3.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA3", "ZDT3");
}

void nsga3_zdt6()
{
    ZDT6 zdt6(10);

    RCGA GA(zdt6.num_vars, zdt6, zdt6.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.8 });
    GA.mutation_method(mutation::real::Gauss{ 1.0 / zdt6.bounds().size() });

    benchmarkMoga(GA, 250, "NSGA3", "ZDT6");
}

void nsga3_dtlz1()
{
    DTLZ1 dtlz1(7, 3);

    RCGA GA(dtlz1.num_vars, dtlz1, dtlz1.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / dtlz1.bounds().size() });

    benchmarkMoga(GA, 1500, "NSGA3", "DTLZ1");
}

void nsga3_dtlz2()
{
    DTLZ2 dtlz2(12, 3);

    RCGA GA(dtlz2.num_vars, dtlz2, dtlz2.bounds());

    GA.population_size(100);
    GA.selection_method(selection::multi_objective::NSGA3{});
    GA.crossover_method(crossover::real::SimulatedBinary{ 0.9, 15.0 });
    GA.mutation_method(mutation::real::Uniform{ 1.0 / dtlz2.bounds().size() });

    benchmarkMoga(GA, 1500, "NSGA3", "DTLZ2");
}

#endif // !NSGA3_BENCHMARK_HPP