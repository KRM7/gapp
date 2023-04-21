/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "genetic_algorithm.hpp"

using namespace genetic_algorithm;
using namespace genetic_algorithm::crossover;
using namespace genetic_algorithm::problems;
using namespace Catch;


TEMPLATE_TEST_CASE("binary_crossover", "[crossover][!mayfail]", binary::SinglePoint, binary::TwoPoint, binary::Uniform)
{
    using Crossover = TestType;

    Rastrigin objective_function(10);
    BinaryGA GA(400);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(Crossover{ 0.75 });
    GA.mutation_method(mutation::binary::Flip{ 0.015 });

    GA.solve(objective_function, 1000);

    double best_found = GA.solutions()[0].fitness[0];
    double maximum = objective_function.optimal_value()[0];

    REQUIRE(best_found <= maximum);
    REQUIRE(best_found == Approx(maximum).margin(1E-6));
}

TEMPLATE_TEST_CASE_SIG("binary_npoint_crossover", "[crossover][!mayfail]", ((size_t N), N), 1, 2, 3, 15)
{
    Rastrigin objective_function(10);
    BinaryGA GA(400);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::binary::NPoint{ 0.75, N });
    GA.mutation_method(mutation::binary::Flip{ 0.015 });

    GA.solve(objective_function, 1000);

    double best_found = GA.solutions()[0].fitness[0];
    double maximum = objective_function.optimal_value()[0];

    REQUIRE(best_found <= maximum);
    REQUIRE(best_found == Approx(maximum).margin(1E-6));
}

TEMPLATE_TEST_CASE("real_crossover", "[crossover][!mayfail]", real::Arithmetic, real::BLXa, real::SimulatedBinary, real::Wright)
{
    using Crossover = TestType;

    Rastrigin objective_function(10);
    RCGA GA(100);

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{} });
    GA.crossover_method(Crossover{ 0.6 });
    GA.mutation_method(mutation::real::NonUniform{ 1.0 / objective_function.num_vars() });

    GA.solve(objective_function, objective_function.bounds(), 1000);

    double best_found = GA.solutions()[0].fitness[0];
    double maximum = objective_function.optimal_value()[0];

    REQUIRE(best_found <= maximum);
    REQUIRE(best_found == Approx(maximum).margin(1E-6));
}

TEMPLATE_TEST_CASE("permutation_crossover", "[crossover][!mayfail]", perm::Order1, perm::Order2, perm::Position, perm::PMX, perm::Cycle, perm::Edge)
{
    using Crossover = TestType;

    TSP52 objective_function;
    PermutationGA GA(500);

    GA.algorithm(algorithm::SingleObjective{ selection::Boltzmann{} });
    GA.crossover_method(Crossover{ 0.9 });
    GA.mutation_method(mutation::perm::Inversion{ 0.95 });

    GA.solve(objective_function, 1000);

    double best_found = GA.solutions()[0].fitness[0];
    double maximum = objective_function.optimal_value()[0];

    REQUIRE(best_found <= maximum);
    REQUIRE(best_found >= 1.2 * maximum);
}

TEMPLATE_TEST_CASE("integer_crossover", "[crossover][!mayfail]", integer::SinglePoint, integer::TwoPoint, integer::Uniform)
{
    using Crossover = TestType;

    StringFinder objective_function("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt.");
    IntegerGA GA(100);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(Crossover{ 0.8 });
    GA.mutation_method(mutation::integer::Uniform{ 0.01 });

    GA.solve(objective_function, objective_function.bounds(), 1000);

    double best_found = GA.solutions()[0].fitness[0];
    double maximum = objective_function.optimal_value()[0];

    REQUIRE(best_found <= maximum);
    REQUIRE(best_found >= 0.95 * maximum);
}

TEMPLATE_TEST_CASE_SIG("integer_npoint_crossover", "[crossover][!mayfail]", ((size_t N), N), 1, 2, 3, 15)
{
    StringFinder objective_function("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Pellentesque gravida ut ipsum at tincidunt.");
    IntegerGA GA(100);

    GA.algorithm(algorithm::SingleObjective{ selection::Roulette{} });
    GA.crossover_method(crossover::integer::NPoint{ 0.8, N });
    GA.mutation_method(mutation::integer::Uniform{ 0.01 });

    GA.solve(objective_function, objective_function.bounds(), 1000);

    double best_found = GA.solutions()[0].fitness[0];
    double maximum = objective_function.optimal_value()[0];

    REQUIRE(best_found <= maximum);
    REQUIRE(best_found >= 0.95 * maximum);
}