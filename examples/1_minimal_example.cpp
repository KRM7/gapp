/* Example showing the minimal usage of a real-encoded, single-objective algorithm. */

#include "gapp.hpp" // include everything
#include <iostream>
#include <cmath>

using namespace gapp;

class SinX : public FitnessFunction<RealGene, 1>
{
    FitnessVector invoke(const Chromosome<RealGene>& x) const override { return { std::sin(x[0]) }; }
};

int main()
{
    Population<RealGene> solutions = RCGA{}.solve(SinX{}, Bounds{ 0.0, 3.14 });

    std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0];
}
