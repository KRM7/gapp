/* Example showing the usage of the GAs with different encodings. */

#include "gapp.hpp"

using namespace gapp;

int main()
{
    // binary encoding
    {
        BinaryGA ga;
        ga.solve(problems::Sphere{ 3 });
    }

    // real encoding
    {
        RCGA ga;
        ga.solve(problems::Sphere{ 3 }, Bounds{ -10.0, 10.0 });
    }

    // permutation encoding
    {
        PermutationGA ga;
        ga.solve(problems::TSP52{});
    }

    // integer encoding
    {
        IntegerGA ga;
        ga.solve(problems::StringFinder{ "Hello" }, Bounds<IntegerGene>{ 'A', 'z' });
    }
}
