/* Example showing the usage of the algorithms in the GAs. */

#include "gapp.hpp"

using namespace gapp;

class MyTournamentSelection : public selection::Selection
{
public:
    const CandidateInfo& selectImpl(const GaInfo&, const PopulationView& pop) const override
    {
        const size_t idx1 = rng::randomIndex(pop);
        const size_t idx2 = rng::randomIndex(pop);

        return (pop[idx1].fitness[0] >= pop[idx2].fitness[0]) ? pop[idx1] : pop[idx2];
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
