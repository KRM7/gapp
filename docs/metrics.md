
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Encodings](encodings.md)  
4. [Algorithms](algorithms.md)  
5. [Genetic operators](genetic-operators.md)  
6. [Stop conditions](stop-conditions.md)  
7. **Metrics**  
8. [Miscellaneous](miscellaneous.md)

------------------------------------------------------------------------------------------------

# Metrics

By default, when you run the GA, the `solve()` function returns a set
of pareto optimal solutions, but you won't have any additional information
about how the population evolved over the run.

In order to get more insight about a run, you can set a number of metrics
that will be tracked throughout the run. Each metric tracks a particular attribute of the population,
and the value of the metric will be recorded in each generation.

## Usage

The tracked metrics must be set before a run, using the `track()`
method of the GAs. You can specify any number of metrics in the arguments
of this function:

```cpp
// Track the minimum and maximum of the population's fitness values in each
// generation, for each objective
GA.track(metrics::FitnessMin{}, metrics::FitnessMax{})
// Run the GA
GA.solve(...);
```

After the run, you can access the metrics using the `get_metric<Metric>()` method.
This will return a reference to the given metric, which stores the values for each generation:

```cpp
// Get the metric tracking the minimal fitness values
const auto& fmin = GA.get_metric<metrics::FitnessMin>();

// Print the lowest fitness value of the first objective in the fourth generation
std::cout << fmin[3][0];
```

## Available metrics

There are a number of metrics already implemented by the library in the
`gapp::metrics` namespace, see:

  - `<metrics/fitness_metrics.hpp>` for fitness related metrics
  - `<metrics/distribution_metrics.hpp>` for metrics tracking information about
    the distribution of the population in the objective space
  - `<metrics/misc_metrics.hpp>` for other metrics

All of the metrics that are implemented in the library work for any objective
function, regardless of the number of objectives, but the distribution metrics
are intended to be used for multi-objective optimization problems.


## Custom metrics

If you want to track something that doesn't have a metric already implemented 
for it, it's possible to implement your own metrics. Metrics must be derived
from the `Monitor` class, and implement the `initialize`, and `update` methods:

```cpp
// The second type parameter of Monitor is the type used to store the
// gathered metrics. This will be used as the type of the data_ field.
class MyMetric : public metrics::Monitor<MyMetric, std::vector<double>>
{
    // (optional) Initialize the metric. Called at the start of a run.
    void initialize(const GaInfo& ga) override { data_.clear(); }

    // Update the metric with a new value from the current generation.
    // Called once in each generation.
    void update(const GaInfo& ga) override
    {
        // You can access the current state of the GA through
        // the ga parameter.
        ...
        data_.push_back(...);
    }
};
```

You can use these custom metrics the same way as any other metric that is implemented
by the library.

------------------------------------------------------------------------------------------------
