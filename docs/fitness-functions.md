
1. [Introduction](introduction.md)  
2. **Fitness functions**  
3. [Encodings](encodings.md)  
4. [Algorithms](algorithms.md)  
5. [Genetic operators](genetic-operators.md)  
6. [Stop conditions](stop-conditions.md)  
7. [Metrics](metrics.md)    
8. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Fitness functions

The `solve` methods of the GAs expect a fitness function
as their first arguments. This fitness function defines
the optimization problem that will be solved by the GA.

## Defining a fitness function

Fitness functions have to be implemented as a class derived
either from `FitnessFunctionBase` or `FitnessFunction`. There
are several things that have to be considered when defining a
fitness function, all of which depend on the particular problem 
we are trying to find the solution for.

### Encoding

The first thing that should be considered is the encoding
type, i.e. how the solutions to the problem should be represented
in the population. There are several options provided by the
library, but user-defined encodings can also be used in the cases
where none of these fit the problem well enough.
The representation is determined by the gene type used in the
fitness function and genetic algorithm classes: the solutions
will be encoded as a vector of the specified gene type.

```cpp
template<typename GeneType>
using Chromosome = std::vector<GeneType>;
```

The options for the gene type provided by the library are:

 - BinaryGene
 - RealGene
 - PermutationGene
 - IntegerGene

The gene type is specified as the type parameter of the
`FitnessFunctionBase` and `FitnessFunction` classes.
This gene type will also determine the GA class that has to be
used to find the optimum the fitness function.

See [encodings.md](encodings.md) for more information regarding the encodings.

### Chromosome length

The length of the chromosomes is another parameter that has to be
considered, and has to be specified for the fitness function.
Generally, the chromosome length is going to be equal to the
number of variables in the problem, but this isn't always true.
For example, if the fitness function uses binary-encoding, a
single variable will likely be represented by multiple binary
genes instead of just a single one.

The chromosome length will also determine the base class that
should be used for defining the fitness function. `FitnessFunctionBase`
can be used for every problem, and the chromosome length is specified
in its constructor. `FitnessFunction` can only be used if the chromosome
length is known at compile-time, as the chromosome length is specified
as a template parameter of the class.

```cpp
template<typename GeneType>
class FitnessFunctionBase;

template<typename GeneType, size_t ChromLen>
class FitnessFunction;
```

### Number of objectives

The number of objectives is another important property of the
fitness functions, but it doesn't have to be specified explicitly
in the definition of the fitness function. Instead, it will be
deduced from the legth of the fitness vectors returned by it.
Note that the fitness function always returns a fitness vector,
even for single-objective problems. In the single-objective case
the length of the fitness vectors will simply be 1.

The library assumes that a fitness function will always return
fitness vectors of the same length, which is equal to the number
of objectives. This can't be changed.

The exact number of objectives is not generally relevant to the
GAs, but whether it is single- or multi-objective will determine
which algorithms can be used in the GAs for the fitness function.
See [algorithms.md](algorithms.md) for more information about the
different algorithms.

### Maximization or minimization

The genetic algorithms in the library will try to find the maximum
of the fitness function, which means that the fitness functions
always have to be implemented for maximization. When trying to find
the minimum of some function, the problem can easily be transformed
into a maximization problem by multiplying the function by -1. This
multiplication is not an issue as all parts of the GAs are implemented
in a way so that they can handle negative fitness values (this includes
operators which might not usually support negative values, such as the
roulette selection operator).

### Example

As a simple example, let's say that we are trying to find the minimum
of the function `f(x) = x^2`. The implementation of the fitness
function could be the following in this case:

```cpp
class XSquare : public FitnessFunction</* GeneType = */ RealGene, /* ChromLen = */ 1>
{
    FitnessVector invoke(const Chromosome<RealGene>& x) const override
    {
        return { -x[0] * x[0] };
    }
};
```

A fitness function for a multi-objective problem would be implemented
the same way, with the only difference being the size of the fitness vector
returned by `invoke`.

## Other fitness function properties

There are some other parameters of the fitness function that can
be specified, but these generally don't have to be changed from
their default values. However, more complex fitness functions might
have to set their values different from these defaults in some cases.

### Dynamic fitness functions

By default, the fitness function is assumed to always return the
same fitness vector for a particular chromosome passed to it as an
argument. This assumption is used to prevent unnecessary fitness
function calls, but this optimization would cause incorrect fitness
vectors to be assigned to some solutions if the assumption is false.

For fitness functions where this assumption would not be true, the value
of the `dynamic` parameter in the constructor of `FitnessFunctionBase`
or `FitnessFunction` has to be set to `true`.

### Variable chromosome lengths

The chromosome length parameter of the fitness function has to be
specified even if variable length chromosomes are used to represent
the solutions. In this case the value of this parameter will only be
used to generate the solutions of the initial population, assuming that
the initial population is generated instead of explicitly specified.

### Example

```cpp
// implementation of a dynamic fitness function
class MyFitnessFunction : public FitnessFunction<RealGene, 1>
{
    MyFitnessFunction() : FitnessFunction(/* dynamic = */ true) {}

    FitnessVector invoke(const Chromosome<RealGene>& x) const override;
};
```

## Other concerns

### Numeric issues

The library in general only assumes that the fitness values returned
by the fitness function are valid numbers (i.e. no `NaN` values will
be returned by it).

Whether infinite fitness values are allowed or not depends on the
selection method used in the GA. If the fitness function can return
fitness vectors that contain infinite values, the selection method
(or the algorithm) will have to be chosen accordingly.

Similarly, whether negative fitness values are allowed depends
on the selection method used, but all selection methods implemented
as part of the library are able to handle negative values without
any issues.

### Thread safety

The candidate solutions in the population are evaluated concurrently
in each generation of a run. As a result of this, the implementation
of the `invoke` method in the derived fitness functions must be thread-safe.

------------------------------------------------------------------------------------------------
