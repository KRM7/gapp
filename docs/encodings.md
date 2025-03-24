
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
  |    `MixedGA`    |    MixedGene    | Complex mixed-encoded  |

The encoding also determines which crossover and mutation methods
can be used in the GA as these genetic operators also depend on the
encoding, and are defined for a particular gene type.

All of the GA classes are defined in the `gapp` namespace.

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

The way the candidate solutions to a problem are going be encoded in the GA is
determined by the gene type. For the simple encodings, the representation of
the solutions will include a single chromosome of the given gene type, which
is effectively a vector of genes.
For mixed encodings, the candidate solutions contain multiple chromosomes,
one for each of the component gene types of the mixed gene type.

```cpp
template<typename GeneType>
using Chromosome = std::vector<GeneType>;
```

The candidates also contain some more information in addition to their
chromosomes, like their fitness vectors and constraint violation vectors,
but these are independent of the gene type. The candidate solutions
themselves are represented by the `Candidate` class.

A population is then made up of several candidates encoded in
this way.

```cpp
template<typename GeneType>
using Population = std::vector<Candidate<GeneType>>;
```


### Mixed encodings

The mixed encoded GA class `MixedGA` does not implement a particular
encoding type, but rather a combination of the other simple encoding
types. Both the `MixedGA` class and the `MixedGene` gene type have an
arbitrary number of template parameters, which specify the component
genes of the mixed encoding. These parameters must be one of the simple
encoding types.

For example, in a GA that uses a mixed encoding of binary and real gene
types, the GA and gene type used would be `MixedGA<BinaryGene, RealGene>`
and `MixedGene<BinaryGene, RealGene>`. `BinaryGene` and `RealGene` are
the component gene types of this mixed encoding.

Each of the component genes must be a valid gene type by itself, and
they must also be unique. Note that the component genes being valid
gene types does not mean that they may not be custom gene types, just
that everything has to be defined for these custom gene types as if
they were used directly as a simple encoding for a GA.

As mentioned above, the solutions contain multiple chromosomes in the
case of mixed encodings, with one chromosome for each of the component
gene types. The crossover and mutation operators are applied separately
for each of these components. See [genetic-operators.md](genetic-operators.md)
for more details about the mixed crossover and mutation operators.


### Variable chromosome lengths

The length of the chromosomes is specified as part of the fitness function.
Normally, this will be a constant value, meaning that all the solutions
will have the same chromosome length throughout a run. However, using a
constant chromosome length is not a requirement as long as all parts of the
GA can handle variable lengths. The parts which must be able to do this are
the:

 - fitness function
 - constraints function (if used)
 - crossover operator
 - mutation operator
 - repair function (if used)

The crossover and mutation operators both provide an
`allow_variable_chrom_length()` method that can be used to check if they
support this or not.

The other operators are implemented by the user in all cases, so it is
up to the user to make sure that their implementations are able to
handle these cases if using variable chromosome lengths is required.

## Custom encodings

It is also possible to use a different, new encoding type by defining a
new GA class. In order to do this, you have to:

 - define the gene type that will be used
 - define a specialization for `GaTraits<GeneType>`
 - specialize `is_bounded_gene<GeneType>` if needed
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
used for the mutation operator. In addition to these, the specialization
must also include a static method `randomChromosome`, which returns a
random chromosome of the given gene type.

```cpp
namespace gapp
{
    template<>
    struct GaTraits<MyGeneType>
    {
        using DefaultCrossover = MyCrossover;
        using DefaultMutation  = MyMutation;

        static Probability defaultMutationRate(size_t chrom_len) { return 1.0 / chrom_len; }

        static Chromosome<MyGeneType> randomChromosome(size_t chrom_len);
        // or, if the gene type is bounded:
        //  static Chromosome<MyGeneType> randomChromosome(size_t chrom_len, BoundsView<MyGeneType> bounds);
    };

}; // namespace gapp
```

Specializing the `is_bounded_gene<T>` struct is only needed if the
gene type used should have it's lower and upper bounds specified for
each gene in the `solve` method. The value should be `true` in this case,
and `false` otherwise (which is the default value used in the primary
template).

```cpp
namespace gapp
{
    // This isn't technically needed, since false is the default value
    template<>
    struct is_bounded_gene<MyGeneType> : std::false_type {};

    // Alternatively, if the custom gene type is bounded:
    //  template<>
    //  struct is_bounded_gene<MyGeneType> : std::true_type {};

} // namespace gapp
```

Note that both the `GaTraits`, and the `is_bounded_gene` template
specializations must be defined in the `gapp` namespace, and before
the actual derived GA class is defined.

The actual GA class should be derived from `GA<T>` using the desired gene
type for the type parameter `T`. The derived class doesn't have to implement
anything, but it may optionally implement an `initialize` method:

```cpp
class MyGA : public GA<MyGeneType> {};
```

In addition to everything above, a crossover and a mutation operator will
also have to be defined for this encoding type, as these operators depend
on the encoding type, and the operators included in the library will not
work for new encodings. See [genetic-operators.md](genetic-operators.md)
for more details about how to define new crossover and mutation methods.

### Mixed gene types

For mixed gene types, the specializations of the `GaTraits` and the `is_bounded_gene`
struct, and also the crossover and mutation methods only have to exist for
the component gene types, not for the actual `MixedGene<>` instance itself.
This means that generally, nothing has to be defined when using some mixed
encoding type, unless it contains a custom gene type as one of its component
gene types, in which case these must be defined for that custom gene type.

The GA class used for the mixed encodings is always an instance of the
`MixedGA` class, even if it contains a custom component gene. For example,
if the component genes of the mixed encoding were `BinaryGene`, `RealGene`,
and `CustomGene`, the GA class would be `MixedGA<BinaryGene, RealGene, CustomGene>`.

------------------------------------------------------------------------------------------------

<p align="right"><a href="algorithms.md">Next: Algorithms</a></p>
