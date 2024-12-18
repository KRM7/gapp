/*
* Example showing a simple real-encoded, constrained optimization problem.
* 
* The objective function is f(x, y) = x^2 + y^2,
* where both x and y are in the closed interval [-1.0, 1.0],
* with 2 additional constraints:
*   c0: x > 0
*   c1: y > 0
* 
* Without the constraints, the function has 4 maximums in the interval
* at (-1.0, -1.0), (-1.0, 1.0), (1.0, -1.0), and (1.0, 1.0).
* Out of these 4, only (1.0, 1.0) satisfies both of the constraints.
*/

#include "gapp.hpp"
#include <iostream>
#include <format>
#include <cmath>

using namespace gapp;

class XYSquare : public FitnessFunction<RealGene, 2>
{
    FitnessVector invoke(const Candidate<RealGene>& sol) const override
    {
        // Compute the value of the objective function f(x, y) = x^2 + y^2
        const double fx = std::pow(sol[0], 2) + std::pow(sol[1], 2);

        // Compute a penalty based on the constraint violation values
        const double cv = sol.constraint_violation[0] + sol.constraint_violation[1];

        return { fx - cv };
    }
};

int main()
{
    RCGA ga{ /* population_size = */ 100 };
    
    /* Specify the constraints. We need to assign some positive constraint violation
       value to negative values of the variables for both constraints. */
    ga.constraints_function([](const GaInfo&, const Chromosome<RealGene>& sol)
    {
        return CVVector{ -sol[0], -sol[1] };
    });

    /* Attempt to improve solutions that violate one of the constraints using a repair function. 
       Note that this part is not necessary for the constraint handling and could just be skipped. */
    ga.repair_function([](const GaInfo&, const Candidate<RealGene>& sol, Chromosome<RealGene>& chrom)
    {
        if (!sol.has_constraint_violation()) return false;

        chrom[0] = -chrom[0];
        chrom[1] = -chrom[1];
        return true;
    });

    /* Run the algorithm and print the results. */
    Population<RealGene> solutions = ga.solve(XYSquare{}, Bounds{ -1.0, 1.0 });

    std::cout << std::format("The maximum of f(x, y) = (x^2 + y^2) in [-1.0, 1.0], "
                             "with the constraints (x > 0) and (y > 0) is at (x = {}, y = {}).\n",
                             solutions[0][0], solutions[0][1]);
}
