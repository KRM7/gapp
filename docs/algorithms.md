
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Constraint handling](constraint-handling.md)  
4. [Encodings](encodings.md)  
5. **Algorithms**  
6. [Genetic operators](genetic-operators.md)  
7. [Stop conditions](stop-conditions.md)  
8. [Metrics](metrics.md)  
9. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Algorithms

The *algorithms* in the library refer to a part of the GAs, which are used
to define different genetic algorithm variants. They consist of a selection and
a population replacement strategy, which define the overall evolution process in
combination with the genetic operators. For some algorithms, the selection and
replacement methods may be specified independently of eachother, while other
algorithms might not allow for this and simply use the methods defined as part
of that algorithm.

The algorithms implemented by the library generally belong to 2 categories:
single-, and multi-objective algorithms. These can only be used for single-
and multi-objective optimization problems respectively. It is also possible
to implement more general algorithms that work for any type of problem, but
the library currently doesn't include such an algorithm.

There are 3 algorithms implemented by the library, all defined in the `gapp::algorithm`
namespace:

 - SingleObjective (single-objective)
 - NSGA-II	(multi-objective)
 - NSGA-III	(multi-objective)

## Selecting the algorithm

By default, if no algorithm is explicitly specified for the GA, a default one will
be automatically selected based on the number of objectives of the fitness function.
As a result of this, the default algorithm used by the GA will always be compatible
with the fitness function regardless of the number of objectives.

```cpp
BinaryGA ga;
ga.solve(f); // run using the default algorithm
```

It is also possible to explicitly select the algorithm to be used by the GA.
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

If an algorithm was explicitly specified, it can be reset by passing a `nullptr`
to the `algorithm` setter. This will result in the default algorithm being used,
the same as in the case where no algorithm was explicitly set.

```cpp
ga.algorithm(nullptr);
ga.solve(f); // uses the default algorithm
```

## The single-objective algorithm

The `SingleObjective` algorithm is not a concrete algorithm implementation
like the NSGA-II and NSGA-III algorithms are. It is simply a wrapper that
combines a selection and a population replacement method, which can be
specified independently of eachother.

The library implements several selection and population replacement methods
that can be used. These are defined in the `gapp::selection` and `gapp::replacement`
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
derived from `algorithm::Algorithm`. The class technically only has to implement
the `selectImpl` and `nextPopulationImpl` methods, but more complex, and efficient
algorithm implementations will want to implement several additional methods.

```cpp
class MyAlgorithm : public algorithm::Algorithm
{
public:
    // ...
};
```

## Custom selection and replacement methods (single-objective)

For the `SingleObjective` algorithms, it is possible to define additional selection
and replacement methods separately without having to define a completely new
algorithm.

Simple selection and population replacement methods can be defined using a lambda
function or some other callable type. As an example, a simple tournament selection
method could be implemented the following way:

```cpp
algorithm::SingleObjective algorithm;

algorithm.selection_method([](const GaInfo& context, const PopulationView& pop)
{
    const size_t idx1 = rng::randomIdx(pop);
    const size_t idx2 = rng::randomIdx(pop);

    return (pop[idx1].fitness[0] >= pop[idx2].fitness[0]) ? idx1 : idx2;
});
```

More complex operators can be implemented as classes derived from `selection::Selection`
and `replacement::Replacement` respectively. The implementation of the simple
tournament selection shown above could also be implemented this way:

```cpp
class MyTournamentSelection : public selection::Selection
{
public:
    size_t selectImpl(const GaInfo& context, const PopulationView& pop) const override
    {
        const size_t idx1 = rng::randomIdx(pop);
        const size_t idx2 = rng::randomIdx(pop);
        
        return (pop[idx1].fitness[0] >= pop[idx2].fitness[0]) ? idx1 : idx2;
    }
};
```

This selection method could then be set for the single-objective algorithm the same
way as in the case of using a lambda function:

```cpp
algorithm::SingleObjective algorithm;
algorithm.selection_method(MyTournamentSelection{});
```

> [!TIP]
> A more general version of the tournament selection operator is
> already implemented in the library by the `selection::Tournament` class.

------------------------------------------------------------------------------------------------

<p align="right"><a href="genetic-operators.md">Next: Genetic operators</a></p>
