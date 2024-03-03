
1. [Introduction](introduction.md)  
2. [Fitness functions](fitness-functions.md)  
3. [Encodings](encodings.md)  
4. [Algorithms](algorithms.md)  
5. [Genetic operators](genetic-operators.md)  
6. [Stop conditions](stop-conditions.md)  
7. [Metrics](metrics.md)  
8. **Miscellaneous**  

------------------------------------------------------------------------------------------------

# Miscellaneous

## Floating-point context

There are multiple places during the runs of the genetic algorithms
where floating-point numbers are compared. First, the fitness vectors
of the solutions need to be compared to determine which solutions
are better. Additionally, when the GA uses real-encoding, the chromosomes
are also encoded as vectors of floating point numbers, so comparing the
candidate solutions can also involve comparing floating point numbers.

These comparisons are not done as exact comparisons, but instead use
an absolute and a relative tolerance value. The actual tolerance used for
a comparison will be the greater of these two tolerances. The values used
for them can be found using the `math::Tolerances::abs()` and
`math::Tolerances::rel()` functions.

```cpp
std::cout << "The default absolute tolerance used is " << math::Tolerances::abs() << "\n";
std::cout << "The default relative tolerance around 1.0 is " << math::Tolerances::rel(1.0) << "\n";
```

The tolerance values can be changed using the `ScopedTolerances`
class, which expects the new absolute and relative tolerance values
in its constructor. These values will be used for the lifetime of
the `ScopedTolerances` instance, and the destructor of the class
will reset the tolerances to their old values.

```cpp
math::ScopedTolerances _(/* abs = */ 1e-10, /* rel = */ 1e-12);
```

Exact comparisons can be used by setting both tolerance values to 0.

```cpp
math::ScopedTolerances _(0.0, 0.0);
```

Note that the tolerance values are global variables which will be used for
the comparisons on every thread, so they should not be modified on multiple
threads concurrently. This means that instances of the `ScopedTolerances`
class should not exist on multiple threads at once.


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

------------------------------------------------------------------------------------------------
