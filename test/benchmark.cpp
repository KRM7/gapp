#include "binary_soga.hpp"
#include "real_soga.hpp"
#include "perm_soga.hpp"
#include "integer_soga.hpp"
#include "nsga2.hpp"
#include "nsga3.hpp"

#include <iostream>
#include "../src/utility/iterators.hpp"

int main()
{
    //binary_rastrigin();
    //binary_rosenbrock();
    //binary_schwefel();
    //binary_griewank();
    //binary_ackley();

    //real_rastrigin();
    //real_rosenbrock();
    //real_schwefel();
    //real_griewank();
    //real_ackley();

    //perm_tsp52();
    //perm_tsp76();
    //perm_tsp124();
    //perm_tsp152();
    //perm_tsp226();
    //perm_tsp299();
    //perm_tsp439();

    //integer_hello();
    //integer_sentence();

    nsga2_kur();
    nsga2_zdt2();
    nsga2_zdt3();
    nsga2_zdt6();

    nsga3_kur();
    nsga3_zdt2();
    nsga3_zdt3();
    nsga3_zdt6();
    
    nsga2_dtlz1();
    nsga2_dtlz2();
    //begin:
    nsga3_dtlz1();
    nsga3_dtlz2();
    //goto begin;
}