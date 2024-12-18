
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Constraint handling](constraint-handling.md)  
4. **Encodings**  
5. [Algorithms](algorithms.md)  
6. [Genetic operators](genetic-operators.md)  
7. [Stop conditions](stop-conditions.md)  
8. [Metrics](metrics.md)  
9. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Encodings

The library defines several GA classes instead of using just a single one.
The difference between these is the encoding type used to represent the
problems and their solutions. Each of the GA classes is for a particular
encoding type, and it can only be used for objective functions using the
same encoding type.

The below table lists each of the GA classes in the library, the encoding
(or gene) type used by them, and the problem (or fitness function) type
they can be used for:

  |     GA class    |  Encoding type  |      Problem type      |
  |:---------------:|:---------------:|:----------------------:|
  |   `BinaryGA`    |   BinaryGene    |    Binary-encoded      |
  |     `RCGA`      |    RealGene     | Floating-point-encoded |
  | `PermutationGA` | PermutationGene |     Combinatorial      |
  |   `IntegerGA`   |   IntegerGene   |    Integer-encoded     |

The encoding also determines which crossover and mutation methods
can be used in the GA as these genetic operators also depend on the
encoding, and are defined for a particular gene type.

All of these GA classes are defined in the main `gapp` namespace.

```cpp
// The fitness function uses permutation encoding, so we
// use PermutationGA
PermutationGA{}.solve(problems::TSP52{});

// The fitness function uses real-encoding, so we use RCGA
RCGA{}.solve(problems::Sphere{}, Bounds{ -10.0, 10.0 });

// The fitness function uses binary-encoding, so we use the
// BinaryGA class
BinaryGA{}.solve(problems::Sphere{});
```


## Solution representation

How the candidate solutions to a problem are going be encoded in the GA is
determined by the gene type. The representation of the solutions will always
be a vector of the gene type used. There is currently no way to change this
to use some other data structure instead of a vector, so this should be taken
into account when defining new encodings.

```cpp
template<typename GeneType>
using Chromosome = std::vector<GeneType>;
```

The candidates contain some more information in addition to their
chromosomes, like their fitness vectors, but these are independent
of the gene type. They are represented by the `Candidate` class.

A population is then made up of several candidates encoded in
this way.

```cpp
template<typename GeneType>
using Population = std::vector<Candidate<GeneType>>;
```


### Mixed encodings

For mixed encodings, where each gene might have a different type, the
recommendation is to use `std::variant` or a similar type for the gene
type.


### Variable chromosome lengths

The length of the chromosomes is specified as part of the fitness function.
Normally, this will be a constant value, meaning that all the solutions
will have the same chromosome length throughout a run. However, using a
constant chromosome length is not a requirement as long as all parts of the
GA can handle variable lengths. The parts which must be able to do this are
the:

 - fitness function
 - constraints function
 - crossover operator
 - mutation operator
 - repair function

The crossover and mutation operators both provide an
`allow_variable_chrom_length()` method that can be used to check if they
support this or not.


## Custom encodings

It is also possible to use a different encoding type by defining a
new GA class. In order to do this, you have to:

 - define the gene type that will be used
 - define a specialization for `GaTraits<GeneType>`
 - specialize `is_bounded<GeneType>` if needed
 - define the GA class, derived from `GA<GeneType>`
 - define crossover and mutation operators for the new encoding

The gene type may be anything, with one restriction: the types
already used for the existing encodings are reserved and can't
be used to define new encodings. See `<encoding/gene_types.hpp>`
for the types that are already in use.

```cpp
using MyGeneType = std::variant<double, int>;
```

The specialization of `GaTraits<T>` for the gene type is required
in order to define some attributes of the GA. These are the default
crossover and mutation operators that will be used by the GA when
they are not specified explicitly, and the default mutation probability
used for the mutation operator:

```cpp
namespace gapp
{
    template<>
    struct GaTraits<MyGeneType>
    {
        using DefaultCrossover = MyCrossover;
        using DefaultMutation  = MyMutation;
        static Probability defaultMutationRate(size_t chrom_len) { return 1.0 / chrom_len; }
    };

}; // namespace gapp
```

Specializing the `is_bounded<T>` variable template is only needed if
the gene type used should have it's lower and upper bounds specified
for each gene in the `solve` method. The value should be `true` in
this case, and `false` otherwise (which is the default value used in
the primary template).

```cpp
namespace gapp
{
    // This isn't technically needed, since false is the default value
    template<>
    inline constexpr bool is_bounded<MyGeneType> = false;

} // namespace gapp
```

The actual GA class should be derived from `GA<T>` using the desired gene
type for the type parameter `T`. The derived class only has to implement
the `generateCandidate` method, and optionally, the `initialize` method:

```cpp
class MyGA : public GA<MyGeneType>
{
public:
    // Generate a random candidate solution. This is used to
    // create the initial population.
    Candidate<GeneType> generateCandidate() const override
    {
        Candidate<GeneType> candidate;
        // ...
        return candidate;
    }
};
```

In addition to everything above, a crossover and a mutation operator will
also have to be defined for this encoding type, as these operators depend
on the encoding type, and the operators included in the library will not
work for new encodings.

------------------------------------------------------------------------------------------------

<p align="right"><a href="algorithms.md">Next: Algorithms</a></p>
