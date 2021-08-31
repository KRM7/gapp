/* 
* Example showing the usage of the GA with a custom gene type.
* Some parts of this example are incomplete, it is not meant to be compiled, but only to show how to create
* a genetic algorithm with a custom encoding type.
*/

#include "../src/base_ga.h" /* For the base GA class. */
#include "../src/rng.h"     /* Some functions for thread-safe random number generation. */

#include <vector>
#include <utility>
#include <functional>
#include <cstddef>

using namespace std;
using namespace genetic_algorithm;


/* Define the gene type that will be used in the algorithm. */
struct GeneType
{
    int first = 0;
    double second = 0.0;
    /* ... */

    /* The equal comparison operator must be defined for the gene type. */
    bool operator==(const GeneType& rhs) const
    {
        return this->first == rhs.first && this->second == rhs.second;
    }
};

/* std::hash must be defined for the gene type used in the GA. */
namespace std
{
    template<>
    struct hash<GeneType>
    {
        size_t operator()(const GeneType& gene) const noexcept
        {
            size_t seed = std::hash<int>()(gene.first);
            seed ^= std::hash<double>()(gene.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
}

/* Implement the genetic algorithm for this encoding. */
class customGA : public GA<GeneType>
{
public:

    /* Constructor. */
    customGA(size_t chrom_len, fitnessFunction_t fitness_function)
        : GA(chrom_len, fitness_function)
    {
    }

private:

    /*
    * Need to override the generateCandidate, crossover, and mutate functions.
    * All of these functions should be thread-safe.
    */

    Candidate generateCandidate() const override
    {
        Candidate sol;

        /* ... */

        return sol;
    }

    CandidatePair crossover(const Candidate& p1, const Candidate& p2) const override
    {
        auto [c1, c2] = tie(p1, p2);

        /* ... */

        return { c1, c2 };
    }

    void mutate(Candidate& child) const override
    {
        /* ... */
    }
};

/* The fitness function used in the algorithm. */
std::vector<double> fitnessFunction(const std::vector<GeneType>& chrom)
{
    double fitness = 0.0;

    /* ... */

    return { fitness };
}

int main()
{
    /* The usage of the GA is the same as the usage of the already implemented ones (eg. BinaryGA, RCGA, etc.) */

    size_t chrom_len = 5;
    customGA GA(chrom_len, fitnessFunction);

    /* ... */

    auto sols = GA.run();

    /* ... */

    return 0;
}