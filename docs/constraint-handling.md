
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. **Constraint handling**  
4. [Encodings](encodings.md)  
5. [Algorithms](algorithms.md)  
6. [Genetic operators](genetic-operators.md)  
7. [Stop conditions](stop-conditions.md)  
8. [Metrics](metrics.md)  
9. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Constrained optimization

Optimization problems frequently have one or more constraints associated
with the objective function. The library provides a straightforward way
to specify an arbitrary number of constraints, but the handling of the
constraints is generally up to the user, due to the many different ways
constraint handling could be implemented.


## Defining constraints

The constraints can be specified using the `constraints_function` method
of the GAs. This constraints function takes a candidate as its parameter,
and returns a vector of constraint violation values for the candidate.

```cpp
using ConstraintsFunction = std::function<CVVector(const GaInfo&, const Candidate<T>&)>;
```

The returned vector should consist of a constraint violation value for each of
the constraints associated with the optimization problem. These values specify
the degree of violation for the constraints. Higher values mean greater degrees
of constraint violation, while 0 or lower values mean that the solution does not
violate that particular constraint.

```cpp
// Specify a constraint that the sum of the variables in the chromosome must be
// greater than 0
ga.constraints_function([](const GaInfo&, const Candidate<RealGene>& sol)
{
    return { -std::accumulate(sol.begin(), sol.end(), 0.0) };
});
```

The constraints are evaluated before the fitness function, and the computed
constraint violation vector is a member of the candidates, so the constraint
violation values are available in the fitness function and also in other parts
of the genetic algorithms.

This means that the constraint violation values can be considered during the
fitness calculation, the selection process, in the repair function call, or any
other part of the algorithms. The only exception is the mutation operator.

The constraint function will be evaluated concurrently for multiple candidate
solutions of the population, so its implementation should be thread-safe.


## Constraint handling methods

There are several ways that the constraint violations can be handled in the
algorithms. Implementing this is up to the user, but we will look at some
examples for some of the possible methods here.

### Constraints as objectives

A simple approach is to consider each constraint as an additional objective,
by including each constraint in the fitness function definition as a separate
dimension. In this case, specifying a constraints function separately is not
necessary, since it's already part of the fitness function definition. The
problem is then solved as an unconstrained multi-objective optimization
problem.

### Penalty based methods

Another method is to consider the constraint violation values in the fitness
function implementation, and penalize the fitness values of the solutions which
violate some constraint. This could be a fixed penalty applied to solutions that
violate a constraint, or a value based on the degree of violation, assigning
larger penalties for greater degrees of violation. Another possibility is to
simply assign the lowest possible fitness value to solutions that violate the
constraints:

```cpp
class Fx : public FitnessFunction</* GeneType = */ RealGene, /* ChromLen = */ 1>
{
    FitnessVector invoke(const Candidate<RealGene>& x) const override
    {
        if (x.has_constraint_violation())
            return { std::numeric_limits<double>::lowest() };

        // return normal fitness ...
    }
};
```

### Repair and other operators

The constraint violations can also be considered in other parts of the GA. For
example, the selection or population replacement methods could be implemented
to prefer solutions that do not violate any constraints. An implementation of
the tournament selection operator which also takes the constraint violations
into account may look like this:

```cpp
class ConstrainedTournamentSelection : public selection::Selection
{
public:
    size_t selectImpl(const GaInfo& context, const PopulationView& pop) const override
    {
        const size_t idx1 = rng::randomIdx(pop);
        const size_t idx2 = rng::randomIdx(pop);

        if (pop[idx1].has_constraint_violation() && !pop[idx2].has_constraint_violation())
            return pop[idx1];
        
        return (pop[idx1].fitness[0] >= pop[idx2].fitness[0]) ? idx1 : idx2;
    }
};
```

It is also possible to specify a repair function that attempts to modify
solutions that violate a constraint so that it does not do so, or at
least in a way that its degree of violation will be smaller:

```cpp
ga.repair_function([](const GaInfo& ga, Candidate<RealGene>& sol)
{
    if (sol.has_constraint_violation())
    {
        // modify the chromosome of the candidate ...
        return true;
    }
    return false;
});
```

These approaches can also be combined with each other, and also with the fitness
penalty based approaches.


## Example

A complete example for solving a constrained optimization problem can be found
[here](../examples/4_constrained_problem.cpp). The example uses a fitness penalty
and a repair function to deal with the constraints.

------------------------------------------------------------------------------------------------

<p align="right"><a href="encodings.md">Next: Encodings</a></p>
