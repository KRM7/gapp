/* Example showing how to change the tolerances used for floating-point comparisons in the GA. */

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
    std::cout << "The default absolute tolerance used is " << math::Tolerances::abs() << "\n";
    std::cout << "The default relative tolerance around 1.0 is " << math::Tolerances::rel(1.0) << "\n\n";

    {
        auto solutions = RCGA{}.solve(SinX{}, Bounds{ 0.0, 3.14 });
        std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0];
        std::cout << ",\twith abs_tol = " << math::Tolerances::abs() << ",\trel_tol = " << math::Tolerances::rel(1.0) << "\n";
    }


    {
        math::ScopedTolerances _(/* abs = */ 0.1, /* rel = */ 0.1);

        auto solutions = RCGA{}.solve(SinX{}, Bounds{ 0.0, 3.14 });
        std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0];
        std::cout << ",\twith abs_tol = " << math::Tolerances::abs() << ",\trel_tol = " << math::Tolerances::rel(1.0) << "\n";
    }


    {
        math::ScopedTolerances _(/* abs = */ 0.0, /* rel = */ 0.0);

        auto solutions = RCGA{}.solve(SinX{}, Bounds{ 0.0, 3.14 });
        std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0];
        std::cout << ",\twith abs_tol = " << math::Tolerances::abs() << ",\trel_tol = " << math::Tolerances::rel(1.0) << "\n";
    }
}
