
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Encodings](encodings.md)  
4. [Algorithms](algorithms.md)  
5. [Genetic operators](genetic-operators.md)  
6. **Stop conditions**  
7. [Metrics](metrics.md)  
8. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Stop conditions

The stop condition of the GAs determine when the run will be terminated.
By default, a run will only stop when reaching the maximum number of generations
specified in `solve()`, or in the case where no maximum was specified there,
it will be the number of generations set using `max_gen()`, or the default
value:

```cpp
PermutationGA GA;
GA.solve(f); // uses the default value of max_generations for stopping

GA.max_gen(1000);
GA.solve(f); // runs the GA for 1000 generations

GA.solve(f, 500); // runs the GA for 500 generations
```

## Early stopping

In a lot of cases, you might want to use a different stopping criterion instead
of using a set number of generations. For example, it could make sense to stop
the run once the solutions have stopped improving significantly, in order to save
time.

It is possible to customize the stop condition of the GAs by specifying
an early-stop condition in addition to the maximum number of generations.

The early-stop condition can be used to terminate a run before it reaches the
maximum number of generations set, when the early-stop condition is met. Note
that even if such a stop condition is set, the GA will still respect the
maximum number of generations set, and it will always stop when reaching that,
regardless of the early-stop condition.

## Usage

There are a number of generally useful stop conditions already implemented
by the library in the `gapp::stopping` namespace. See
`<stop_condition/stop_condition.hpp>` for more details.

The early-stop condition can be set using the `stop_condition()` method:

```cpp
BinaryGA GA;
// stop the run once the average fitness values of the population
// haven't improved significantly for 5 generations
GA.stop_condition(stopping::FitnessMeanStall{ 5 });
GA.solve(f);
```

If you want to clear the early-stop condition, you can either set it to a
`nullptr`, or use `stopping::NoEarlyStop`:

```cpp
GA.stop_condition(nullptr);
// or
GA.stop_condition(stopping::NoEarlyStop{});
```

## Composite early-stop conditions

Early-stop conditions can be combined to create more complex stop conditions
using `operator&&`, `operator||`, and `operator!`. This allows for specifying more complex
stopping criterion without having to write your own custom stop conditions:

```cpp
GA.stop_condition(stopping::FitnessMeanStall{ 5 } && stopping::FitnessBestStall{ 5 });
```

## Custom early-stop conditions

If the stop conditions already implemented by the library are not enough,
you can also define your own stop conditions.

For simple stop conditions, you can use a lambda function (or any other callable):

```cpp
GA.stop_condition([](const GaInfo& ga)
{
    // Check if the stopping criterion is met,
    // return true if the run should be terminated
    ...
    return false;
});
```

For more complex stop conditions, you can define your own stop condition class.
The class should be derived from `stopping::StopCondition`, and implement
`stop_condition`, and optionally `initialize`:

```cpp
class MyStopCondition : public stopping::StopCondition
{
    // Check if the stopping criterion is met.
    // Called once in every generation
    bool stop_condition(const GaInfo& ga) override
    {
        ...
        // Return true if the run should be terminated
        return false;
    }

    // Initialize the stop condition. Called once at the start
    // of a run. Implementing this is optional.
    void initialize(const GaInfo& ga) override { ... }
};
```

------------------------------------------------------------------------------------------------
