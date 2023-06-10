/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "binary_soga.hpp"
#include "real_soga.hpp"
#include "perm_soga.hpp"
#include "integer_soga.hpp"
#include "nsga2.hpp"
#include "nsga3.hpp"

int main()
{
    //while (true)
    {
        binary_sphere();
        binary_rastrigin();
        binary_rosenbrock();
        binary_schwefel();
        binary_griewank();
        binary_ackley();
        binary_levy();

        real_sphere();
        real_rastrigin();
        real_rosenbrock();
        real_schwefel();
        real_griewank();
        real_ackley();
        real_levy();

        perm_tsp52();
        perm_tsp76();
        perm_tsp124();
        perm_tsp152();
        perm_tsp226();
        perm_tsp299();
        perm_tsp439();

        integer_hello();
        integer_sentence();

        benchmark_nsga2_zdt();
        benchmark_nsga2_dtlz();

        benchmark_nsga3_zdt();
        benchmark_nsga3_dtlz();
    }
}