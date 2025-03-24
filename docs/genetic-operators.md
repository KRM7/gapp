
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Constraint handling](constraint-handling.md)  
4. [Encodings](encodings.md)  
5. [Algorithms](algorithms.md)  
6. **Genetic operators**  
7. [Stop conditions](stop-conditions.md)  
8. [Metrics](metrics.md)  
9. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Genetic operators

The genetic operators are used to create new candidate solutions from the existing
ones in the population, providing the basic search mechanism of the genetic
algorithms. The 2 main operators are the crossover and mutation, but the library
also allows for specifying a repair function.

The selection method is considered to be a part of the *algorithms* in the library,
so they are not discussed here. See [algorithms.md](algorithms.md) for more details
about the selection methods.

The library contains implementations of several crossover and mutation methods
that can be used. These can be found in the `gapp::crossover` and
`gapp::mutation` namespaces respectively.

As the genetic operators operate on candidate solutions, their implementations depend
on the encoding type used for the GA. A particular crossover or mutation method can
therefore only be used with the encoding type it is implemented for.
Because of this, the implemented crossover and mutation methods are further broken
up into multiple namespaces based on the encoding type they can be used with.
For example, the crossover operators are in the following 4 namespaces:

 - `crossover::binary`
 - `crossover::real`
 - `crossover::perm`
 - `crossover::integer`

Crossover methods in the `binary` namespace can only be used for the `BinaryGA`,
methods in the `real` namespace can only be used for the `RCGA`, and so on.
The mutation methods are organized the same way.

The library doesn't provide any repair functions since their use in the GAs
is optional, and their implementation is highly dependent on the problem itself.
These have to be defined by the user if they are needed.

## Crossover

The crossover operator is responsible for generating new solutions from existing
ones. The operator takes 2 solutions selected from the population, and performs
some operation on them with a given probability to generate 2 new solutions.
When the operation is not performed, it returns copies of the parent solutions.

The crossover operator used by the GA can be set either in the constructor or by
using the `crossover_method` method:

```cpp
PermutationGA GA;
GA.crossover_method(crossover::perm::Edge{});
```

The probability of performing the crossover operation is a general parameter
of the crossovers, and it can be set for all of the crossover operators in
their constructors or by using the `crossover_rate` method:

```cpp
PermutationGA GA;
GA.crossover_method(crossover::perm::Edge{ /* crossover_rate = */ 0.8 });
```

The crossover method used by the GA classes can also be accessed using the
`crossover_method` function. This can be used, for example, to change the
crossover probability after the crossover method has already been set:

```cpp
PermutationGA GA;
GA.crossover_method(crossover::perm::Edge{});
GA.crossover_method().crossover_rate(0.8);
```

Note that the parameters of the crossover operators, such as the crossover
probability, are associated with the crossover operators, not the GA itself,
so they will be changed when a new crossover operator is set for the GA.

The crossover operators may also have additional parameters that are specific
to that particular operator.

### Mixed crossovers

The crossover operators for the mixed encodings are composite operators made
up of an independent crossover operator for each of the component gene types
of the mixed encoding. These component operators are applied to each of the
associated chromosomes of the candidates separately.

The crossover class used for all of the mixed crossovers is `crossover::Mixed`.
The component crossovers are specified in the constructor of the class:

```cpp
MixedGA<RealGene, BinaryGene> ga;
ga.crossover_method(crossover::Mixed{ crossover::real::Wright{}, crossover::binary::Uniform{} });
```

The parameters of the crossovers, such as the crossover probability, are associated
with the individual components of the mixed crossover, not with the mixed crossover
itself. This means that each component crossover has its own crossover probability.

The component crossovers can be accessed through the mixed crossover, using the
`component` method. This can be used, for example, to change the crossover
probabilities of the individual components:

```cpp
ga.crossover_method().component<RealGene>().crossover_rate(0.8);
```

There are also some helper methods to access common parameters like the crossover
probability directly:

```cpp
ga.crossover_method().crossover_rate<RealGene>(0.8);
```

## Mutation

The mutation operator is applied to each of the solutions generated by the
crossovers in order to promote diversity in the population. This helps the GA
with exploring more of the search space and avoiding convergence to local
optima.

The mutation operator used by the GAs can be chosen similar to the crossover
operators, either in the constructor of the GA or by using the `mutation_method`
method:

```cpp
RCGA GA;
GA.mutation_method(mutation::real::NonUniform{});
```

Similar to the crossovers, the mutation operators also have a mutation
probability parameter, but how this probability is interpreted depends
on the specific mutation operator. It may either be interpreted on a
per-gene or on a per-solution basis.

The mutation probability can be set similar to the crossover probability,
in the constructor of the mutation operator, or by using the `mutation_rate`
method:

```cpp
RCGA GA;
GA.mutation_method(mutation::real::NonUniform{ /* mutation_rate = */ 0.1 });
```

The mutation method used by the GA class can also be accessed using the
`mutation_method` function. This can be used, for example, to change the
mutation probability after the mutation method has already been set:

```cpp
RCGA GA;
GA.mutation_method(mutation::real::NonUniform{});
GA.mutation_method().mutation_rate(0.1);
```

Note that just like for the crossovers, the mutation rate is associated with the
mutation operator, not the GA itself, which means that it will be changed along
with the mutation operator when that is changed.

Note that just like in the case of the crossover operators, the parameters of
the mutation operators, such as the mutation probability, are associated with
the mutation operators themselves, not with the GA, so they will be changed
when a new mutation operator is set for the GA.

The mutation operators may also have additional parameters that are specific
to that particular operator.

### Mixed mutations

The mixed mutation operators are structured the same as the mixed crossover
operators. They are composite operators made up of an independent mutation
operator for each of the component gene types of the mixed encoding. These
component operators are applied to each of the associated chromosomes of the
candidates separately.

The mutation class used for all of the mixed mutations is `mutation::Mixed`.
The component mutations are specified in the constructor of the class:

```cpp
MixedGA<RealGene, BinaryGene> ga;
ga.mutation_method(mutation::Mixed{ mutation::real::Uniform{}, mutation::binary::Flip{} });
```

The parameters of the mutations, such as the mutation probability, are associated
with the individual components of the mixed mutation, not with the mixed mutation
itself. This means that each component mutation has its own mutation probability.

The component mutations can be accessed through the mixed mutation, using the
`component` method. This can be used, for example, to change the mutation
probabilities of the individual components:

```cpp
ga.mutation_method().component<RealGene>().mutation_rate(0.8);
```

There are also some helper methods to access common parameters like the mutation
probability directly:

```cpp
ga.mutation_method().mutation_rate<RealGene>(0.8);
```

## Repair

The repair function is an additional operator that will be applied to each
solution after the mutations. Using a repair function is optional, and they
are not used in the GAs by default.

The repair function can be specified as some callable object using the `repair_function`
method of the GAs:

```cpp
ga.repair_function([](const GaInfo&, Candidate<RealGene>& sol)
{
    // change the chromosome of the candidate in some way ...

    return true; // return true if the chromosome was changed, or false otherwise
});
```

The implementation of the repair function may not change any part of the
candidates other than their chromosomes. The function should also be thread-safe.

The repair function used by the GA can be cleared by setting it to a nullptr:

```cpp
GA.repair_function(nullptr);
```

## Custom genetic operators

In addition to the operators already implemented in the library,
user defined crossover and mutation operators can also be used in the GAs.

The simplest way to do this is to use a lambda function (or some other callable):

```cpp
RCGA ga;

ga.crossover_method([](const GaInfo&, const Candidate<RealGene>& parent1, const Candidate<RealGene>& parent2)
{
    auto child1 = parent1;
    auto child2 = parent2;

    // perform the crossover ...

    return CandidatePair<RealGene>{ std::move(child1), std::move(child2) };
});

ga.mutation_method([](const GaInfo& ga, const Candidate<RealGene>& sol, Chromosome<RealGene>& chrom)
{
    for (RealGene& gene : chrom)
    {
        if (rng::randomReal() < ga.mutation_rate())
        {
            // modify the gene ...
        }
    }
});
```

Alternatively, crossover and mutation operators can also be implemented as
classes derived from `crossover::Crossover<GeneType>` and
`mutation::Mutation<GeneType>` respectively. Crossovers must implement the
`crossover` method, while mutations must implement the `mutate` method:

```cpp
class MyCrossover : public crossover::Crossover<RealGene>
{
public:
    using Crossover::Crossover;

    CandidatePair<RealGene> crossover(const GaInfo& ga, const Candidate<RealGene>& parent1, const Candidate<RealGene>& parent2) const override
    {
        // perform the crossover ...
    }
};
```

```cpp
class MyMutation : public mutation::Mutation<RealGene>
{
public:
    using Mutation::Mutation;

    void mutate(const GaInfo& ga, const Candidate<RealGene>& candidate, Chromosome<RealGene>& chromosome) const override
    {
        // perform the mutation on chromosome ...
    }
};
```

There are a couple of things that should be kept in mind when implementing these operators
regardless of how they are defined:

 - The crossover implementation should just perform the crossover operation without
   taking the crossover rate into account. Handling the crossover rate is done elsewhere.
 - The mutation implementation must take the mutation rate into account, as the
   interpretation of the mutation rate depends on the specific mutation method
   (it can either be a per-gene or per-chromosome probability).
 - The mutation modifies the `chromosome` parameter directly, and does not return anything.
 - The implementations should be thread-safe. It can be assumed that the mutation operator
   has exclusive access to its `chromosome` parameter.
 - If a mixed encoding uses a custom gene type as one of its component gene types, you
   only need to define a crossover or a mutation operator for that particular custom gene
   type, which can then be used in the mixed crossover or mutation operators.

------------------------------------------------------------------------------------------------

<p align="right"><a href="stop-conditions.md">Next: Stop conditions</a></p>
