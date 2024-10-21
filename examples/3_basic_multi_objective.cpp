/* Example showing a more detailed usage of a real-encoded, multi-objective algorithm. */

#include "gapp.hpp" // include everything
#include <iostream>
#include <cmath>

using namespace gapp;

class Kursawe2 : public FitnessFunction<RealGene, 2>
{
    FitnessVector invoke(const Candidate<RealGene>& sol) const override
    {
        const double f1 = 10.0 * std::exp(-0.2 * std::sqrt(std::pow(sol[0], 2) + std::pow(sol[1], 2)));
        const double f21 = std::pow(std::abs(sol[0]), 0.8) + 5.0 * std::sin(std::pow(sol[0], 3));
        const double f22 = std::pow(std::abs(sol[1]), 0.8) + 5.0 * std::sin(std::pow(sol[1], 3));

        return { -f1, -f21 - f22 };
    }
};

int main()
{
    RCGA GA{ /* population_size = */ 20 };

    GA.algorithm(algorithm::NSGA3{});
    GA.crossover_method(crossover::real::BLXa{});
    GA.mutation_method(mutation::real::Boundary{ 0.1 });
    GA.stop_condition(stopping::FitnessMeanStall{ 5 });

    auto solutions = GA.solve(Kursawe2{}, Bounds{ 0.0, 3.14 });

    std::cout << "The optimal solutions of the Kursawe function in 2 dimensions are:\n x1\t|\tx2\n";
    for (const Candidate<RealGene>& sol : solutions)
    {
        std::cout << sol.chromosome[0] << "\t|\t" << sol.chromosome[1] << "\n";
    }
}
