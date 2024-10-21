
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Constraint handling](constraint-handling.md)  
4. [Encodings](encodings.md)  
5. [Algorithms](algorithms.md)  
6. [Genetic operators](genetic-operators.md)  
7. [Stop conditions](stop-conditions.md)  
8. [Metrics](metrics.md)  
9. **Miscellaneous**  

------------------------------------------------------------------------------------------------

# Miscellaneous

## Floating-point context

There are multiple places during the runs of the genetic algorithms
where floating-point numbers are compared. First, the fitness vectors
of the solutions need to be compared with eachother to determine which solutions
are better. Additionally, when the GA uses real-encoding, the chromosomes
are also encoded as vectors of floating point numbers, so comparing the
candidate solutions will also involve comparing floating point numbers in this
case.

These comparisons are not done as exact comparisons, but instead use
an absolute and a relative tolerance value. The actual tolerance used for
a comparison will be the greater of these two tolerances. The values used
for them can be found using the `math::Tolerances::abs` and
`math::Tolerances::rel` functions.

```cpp
std::cout << "The default absolute tolerance used is " << math::Tolerances::abs() << "\n";
std::cout << "The default relative tolerance around 1.0 is " << math::Tolerances::rel(1.0) << "\n";
```

The tolerance values can be changed using the `ScopedTolerances` class,
which expects the new absolute and relative tolerance values to be specified
in its constructor. These values will be used for the lifetime of the
`ScopedTolerances` instance, and the destructor of the class will reset the
tolerances to their old values.

```cpp
math::ScopedTolerances _(/* abs = */ 1e-10, /* rel = */ 1e-12);
```

Exact comparisons can be used by setting both tolerance values to 0.

```cpp
math::ScopedTolerances _(0.0, 0.0);
```

Note that these tolerance values are global values which will be used for
the comparisons on every thread, and changing their values is not thread-safe.
This means that instances of the `ScopedTolerances` class should not be created
while a genetic algorithm is running, and instances of this class should also
not exist on multiple threads at once.


## Random number generation

Several parts of the GAs depend on random numbers in their
implementations. These numbers are generated using a single
global pseudo-random number generator instance. This PRNG
instance can be accessed as `rng::prng`. There are also
several utility functions for generating random numbers using
this engine in the `rng` namespace, so this generator doesn't
have to be used directly.

The global generator is seeded using a constant value determined
by the value of the `GAPP_SEED` macro. The value of this
can be changed by defining this macro on the command line
while building and using the library.

Alternatively, the PRNG can also be reseeded using its `seed` method.

```cpp
rng::prng.seed(new_seed);
```

The methods of the PRNG and all of the random generation
utilities are thread-safe, and can be used freely by the user
if needed, for example in the implementation of custom
genetic operators. The only exception to this is the `seed()` method
of the prng, which is not thread safe and shouldn't be called
concurrently with the random number generation methods. In practice,
this means that `seed()` sholdn't be called while a GA is running.


## Execution

By default, the library will use multiple threads for running the
genetic algorithms, with the number of threads being the number of
hardware threads available as indicated by `std::thread::hardware_concurrency`.

This can be changed using the `execution_threads` function. The number of threads
used will be the value specified as the argument to this function.

```cpp
execution_threads(1); // run everything on a single thread
```

The specified thread count should be between 1 and the number of hardware threads.
Numbers larger than the number of hardware threads can be set, but they will likely
lead to worse performance. If 0 is specified, it will be ignored and a single thread
will be used instead.

Note that this function is not thread-safe, and shouldn't be called while a genetic
algorithm is running.


## Determinism and reproducibility

Due to the library's reliance on random numbers generated from a global generator, 
the solutions to a problem will generally vary between runs even if the same parameters
are used. However, it is possible to reproduce the results of previous runs by seeding
the random number generator appropriately.

```cpp
rng::prng.seed(0x9e3779b97f4a7c15);
const auto solutions1 = ga.solve(f);

rng::prng.seed(0x9e3779b97f4a7c15);
const auto solutions2 = ga.solve(f);

assert(solutions1 == solutions2);
```

As long as the PRNG is seeded with the same value before each run, the results from the
runs will be the same. This will be true regardless of the number of threads used, there
is no need to use single-threaded execution for this. However, the number of threads used
should not be changed between the runs, as that would lead to different results.

The results also depend on the implementation of the standard library, so they will
not be reproducible using different implementations. This also means that results
are not generally reproducible across different platforms.

------------------------------------------------------------------------------------------------
