
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Encodings](encodings.md)  
4. **Algorithms**  
5. [Genetic operators](genetic-operators.md)  
6. [Stop conditions](stop-conditions.md)  
7. [Metrics](metrics.md)    
8. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Algorithms

The algorithms are used in the library to define different genetic algorithm
variants. They consist of the selection and population replacement methods,
which define the overall evolution process in combination with the genetic
operators.

The algorithms in the library belong to 2 categories: single-, and
multi-objective algorithms. These can only be used for single- and
multi-objective optimization problems respectively. It is possible
to implement general algorithms that work for any type of problem, but the
library currently doesn't have such an algorithm.

There are 3 algorithms provided by the library:

 - SingleObjective (single-objective)
 - NSGA-II	(multi-objective)
 - NSGA-III	(multi-objective)

 All of these algorithms are in the `gapp::algorithm` namespace.

## Selecting the algorithm

By default, if no algorithm is specified for the GA, one will automatically
be selected based on the number of objectives of the fitness function being used.
This means that the default algorithm used by the GA will always be compatible
with the fitness functions regardless of the number of objectives.

```cpp
BinaryGA ga;
ga.solve(f); // run using the default algorithm
```

It is also possible to select a different algorithm to be used by the GA.
This can be done either in the constructor or by using the `algorithm` method:

```cpp
BinaryGA ga;
ga.algorithm(algorithm::NSGA3{});
ga.solve(f);
```

The only thing that should be kept in mind when selecting the algorithm this
way is that it needs to be compatible with the fitness function used.
The single-objective algorithms can only be used for single-objective fitness
functions, and the multi-objective algorithms can only be used with multi-objective
fitness functions.

If an algorithm was explicitly specified, it can be cleared by passing a `nullptr`
to the `algorithm` setter. This will result in the default algorithm being used,
as in the case where no algorithm was explicitly set.

```cpp
ga.algorithm(nullptr);
ga.solve(f); // uses the default algorithm
```

## The single-objective algorithm

The `SingleObjective` algorithm is not a concrete algorithm implementation
like the NSGA-II and NSGA-III algorithms are. It is simply a wrapper that
combines a selection and a population replacement method. These methods
can be selected independently of eachother in the `SingleObjective` algorithm.

The library implements several selection and population replacement methods
that can be used. These are in the `gapp::selection` and `gapp::replacement`
namespaces respectively.

```cpp
BinaryGA ga;
ga.algorithm(algorithm::SingleObjective{}); // use the default selection and replacement methods
ga.solve(f);

// use tournament selection, and elitism for the population replacement methods
ga.algorithm(algorithm::SingleObjective{ selection::Tournament{}, replacement::Elitism{ 5 } });
ga.solve(f);
```

## Custom algorithms

In addition to the algorithms provided by the library, it is also possible to
use user-defined algorithms in the GAs. These must be implemented as a class
that is derived from `algorithm::Algorithm`. The class technically only has to implement
the `selectImpl` and `nextPopulationImpl` methods, but more complex, and efficient
algorithm implementations will have to implement several additional methods.

```cpp
class MyAlgorithm : public algorithm::Algorithm
{
public:
    // ...
};
```

## Custom selection and replacement methods (single-objective)

For the `SingleObjective` algorithms, it's possible to define additional selection
and replacement methods separately without having to define a completely new
algorithm.

Simple selection and population replacement methods can be defined using a lambda
or some other callable type. As an example, a simple tournament selection method
could be implemented this way:

```cpp
algorithm::SingleObjective algorithm;
algorithm.selection_method([](const GaInfo& context, const FitnessMatrix& fmat)
{
    size_t first  = rng::randomIdx(fmat);
    size_t second = rng::randomIdx(fmat);

    return (fmat[first][0] >= fmat[second][0]) ? first : second;
});
```

More complex operators can be implemented as classes derived from `selection::Selection`
and `replacement::Replacement` respectively. The implementation of the simple
tournament selection shown above could also be implemented this way:

```cpp
class MyTournamentSelection : public selection::Selection
{
public:
    size_t selectImpl(const GaInfo& context, const FitnessMatrix& fmat) const override
    {
        size_t first  = rng::randomIdx(fmat);
        size_t second = rng::randomIdx(fmat);
        
        return (fmat[first][0] >= fmat[second][0]) ? first : second;
    }
};
```

Note that a more general version of the tournament selection operator
is already implemented in the library, called `selection::Tournament`.

------------------------------------------------------------------------------------------------
