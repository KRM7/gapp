/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_HPP
#define GA_PROBLEMS_HPP

#include "benchmark_function.hpp"
#include "single_objective.hpp"
#include "multi_objective.hpp"
#include "many_objective.hpp"
#include "travelling_salesman.hpp"
#include "integer.hpp"

/**
* Implementations of some benchmark functions that can be used to test the
* genetic algorithms. There are some benchmark problems implemented for every encoding
* type, and the real-encoded problems can also be used for the binary-encoded %GAs.
* 
* All of the problems are implemented for maximization.
*/
namespace genetic_algorithm::problems {}

#endif // !GA_PROBLEMS_HPP