/* Example showing the usage of the algorithms in the GAs. */

#include "gapp.hpp"

using namespace gapp;

class MyTournamentSelection : public selection::Selection
{
public:
    size_t selectImpl(const GaInfo&, const FitnessMatrix& fmat) const override
    {
        size_t first = rng::randomIndex(fmat);
        size_t second = rng::randomIndex(fmat);

        return (fmat[first][0] >= fmat[second][0]) ? first : second;
    }
};

int main()
{
    BinaryGA ga;
    ga.solve(problems::Sphere{ 3 }); // run using the default algorithm
    ga.solve(problems::Kursawe{});   // the default algorithm works with both single- and multi-objective problems

    // using a different algorithm

    ga.algorithm(algorithm::NSGA3{});
    ga.solve(problems::Kursawe{}); // the NSGA3 algorithm only works with multi-objective problems

    // going back to the default algorithm

    ga.algorithm(nullptr);
    ga.solve(problems::Sphere{ 3 });
    ga.solve(problems::Kursawe{});

    // selecting the selection and replacement methods for the SingleObjective algorithm

    ga.algorithm(algorithm::SingleObjective{ selection::Tournament{}, replacement::Elitism{ 5 } });
    ga.solve(problems::Sphere{ 3 });

    // user defined selection methods

    ga.algorithm(algorithm::SingleObjective{ MyTournamentSelection{} });
    ga.solve(problems::Sphere{ 3 });
}
