/* Example showing the definition of fitness function for the GAs. */

#include "gapp.hpp"
#include <iostream>

using namespace gapp;

// single-objective fitness function
class XSquare : public FitnessFunction</* GeneType = */ RealGene, /* ChromLen = */ 1>
{
public:
    FitnessVector invoke(const Chromosome<RealGene>& x) const override
    {
        return { -x[0] * x[0] };
    }
};

// multi-objective fitness function
class XSquareMulti : public FitnessFunction<RealGene, 1>
{
public:
    FitnessVector invoke(const Chromosome<RealGene>& x) const override
    {
        const double f1 = -x[0] * x[0];
        const double f2 = (std::abs(x[0]) <= 2.0) ? 0.0 : 1.0;
        
        return { f1, f2 };
    }
};

// dynamic fitness function
class XSquareDynamic : public FitnessFunction< RealGene, 1>
{
public:
    XSquareDynamic() : FitnessFunction(/* dynamic = */ true) {}

    FitnessVector invoke(const Chromosome<RealGene>& x) const override
    {
        return { -x[0] * x[0] + rng::randomNormal(0.0, 1.0) };
    }
};

int main()
{
    RCGA ga;

    // single-objective fitness function
    {
        auto solutions = ga.solve(XSquare{}, Bounds{ -100.0, 100.0 });
        std::cout << "The minimum of x^2 in [-100.0, 100.0] is at x = " << solutions[0].chromosome[0] << "\n";
    }

    // multi-objective fitness function
    {
        auto solutions = ga.solve(XSquareMulti{}, Bounds{ -100.0, 100.0 });
        std::cout << "The minimum of x^2 in [-100.0, 100.0] is at x = " << solutions[0].chromosome[0] << "\n";
    }

    // dynamic fitness function
    {
        auto solutions = ga.solve(XSquareDynamic{}, Bounds{ -100.0, 100.0 });
        std::cout << "The minimum of x^2 in [-100.0, 100.0] is at x = " << solutions[0].chromosome[0] << "\n";
    }
}
