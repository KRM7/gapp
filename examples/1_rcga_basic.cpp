/* Example showing the minimal usage of a real-encoded, single-objective algorithm. */

#include "genetic_algorithm.hpp" // include everything
#include <iostream>
#include <cmath>

using namespace genetic_algorithm;

class SinX final : public SingleObjFitnessFunction<RealGene, 1>
{
    FitnessVector invoke(const Chromosome<RealGene>& chrom) const override { return { std::sin(chrom[0]) }; }
};

int main()
{
    SinX objective_function;
    GeneBounds bounds{ 0.0, 3.14 };

    RCGA GA(objective_function, bounds);

    Population<RealGene> solutions = GA.run();

    std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0] << ".\n";
}