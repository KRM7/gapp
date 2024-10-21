/* Example showing a more detailed usage of a real-encoded, single-objective algorithm. */

#include "gapp.hpp" // include everything
#include <iostream>
#include <cmath>

using namespace gapp;

class SinX : public FitnessFunction<RealGene, 1>
{
    FitnessVector invoke(const Candidate<RealGene>& sol) const override { return { std::sin(sol[0]) }; }
};

int main()
{
    RCGA GA{ /* population_size = */ 100 };

    GA.algorithm(algorithm::SingleObjective{ selection::Tournament{}, replacement::KeepBest{} });
    GA.crossover_method(crossover::real::Wright{ 0.8_p });
    GA.mutation_method(mutation::real::Gauss{ 0.1_p });
    GA.stop_condition(stopping::FitnessBestStall{ 5 });

    auto solutions = GA.solve(SinX{}, Bounds{ 0.0, 3.14 });

    std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0];
}
