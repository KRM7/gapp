
# Genetic Algorithms in C++ 

[![linux](https://github.com/KRM7/genetic-algorithms/actions/workflows/linux.yml/badge.svg?branch=master)](https://github.com/KRM7/genetic-algorithms/actions/workflows/linux.yml)
[![windows](https://github.com/KRM7/genetic-algorithms/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/KRM7/genetic-algorithms/actions/workflows/windows.yml)
[![macos](https://github.com/KRM7/gapp/actions/workflows/macos.yml/badge.svg?branch=master)](https://github.com/KRM7/gapp/actions/workflows/macos.yml)
[![sanitizers](https://github.com/KRM7/genetic-algorithms/actions/workflows/sanitizers.yml/badge.svg?branch=master)](https://github.com/KRM7/genetic-algorithms/actions/workflows/sanitizers.yml)
[![code analysis](https://github.com/KRM7/genetic-algorithms/actions/workflows/analysis.yml/badge.svg?branch=master)](https://github.com/KRM7/genetic-algorithms/actions/workflows/analysis.yml)
[![docs](https://github.com/KRM7/gapp/actions/workflows/docs.yml/badge.svg?branch=master)](https://github.com/KRM7/gapp/actions/workflows/docs.yml)

## Overview

gapp is a library of genetic algorithm implementations in C++ for solving single-
and multi-objective optimization problems. The algorithms are highly customizable,
with all of their parts possibly defined by the user, but the library also includes
GAs for several commonly used encoding types, frequently used crossover and mutation
methods for each of these encodings, several stop conditions, and other utilities that
can be used.


## Usage example

```cpp
#include <gapp/gapp.hpp>
#include <iostream>

using namespace gapp;

class SinX : public FitnessFunction<RealGene, 1> 
{
    FitnessVector invoke(const Candidate<RealGene>& x) const override { return { std::sin(x[0]) }; }
};

int main()
{
    auto solutions = RCGA{}.solve(SinX{}, Bounds{ 0.0, 3.14 });

    std::cout << "The maximum of sin(x) in [0.0, 3.14] is at x = " << solutions[0].chromosome[0];
}
```

Possible console output:

```text
The maximum of sin(x) in [0.0, 3.14] is at x = 1.57079
```


## Requirements

The following are needed for building and using the library:

- C++20 compiler (gcc 11, clang 15, msvc 14.30 or later)
- CMake 3.21 or later
- Catch2 3.3 or later (optional, only needed for the tests)


## Installing the library

*See [install-guide.md](docs/install-guide.md) for a more detailed installation and usage guide.*

The library uses CMake as its build system. Assuming you have all of the requirements
listed above, the steps for installing the library in Release config are:

```shell
git clone https://github.com/KRM7/gapp.git --branch v1.0.0
cd gapp/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build . --config Release
sudo cmake --install . --config Release
```

Alternatively, you can also use the install script that is provided with the library, which
will install all available configurations:

```shell
git clone https://github.com/KRM7/gapp.git --branch v1.0.0
sudo bash ./gapp/tools/install.sh
```

Once the library is installed, you can import it using `find_package` and then link
against the `gapp::gapp` target provided by the library:

```cmake
find_package(gapp REQUIRED)
target_link_libraries(YourTarget PRIVATE gapp::gapp)
```


## Documentation

The documentation for the library can be found in the [docs](./docs) directory:

* [Introduction](./docs/introduction.md)
* [Fitness functions](./docs/fitness-functions.md) 
* [Encodings](./docs/encodings.md)  
* [Algorithms](./docs/algorithms.md)  
* [Genetic operators](./docs/genetic-operators.md)  
* [Stop conditions](./docs/stop-conditions.md)  
* [Metrics](./docs/metrics.md)    
* [Miscellaneous](./docs/miscellaneous.md)

The API documentation is available [here](https://krm7.github.io/gapp/).


## Examples

The [examples](./examples) directory contains several examples for using the library.

* [Minimal example](./examples/1_minimal_example.cpp)
* [Single-objective optimization](./examples/2_basic_single_objective.cpp)
* [Multi-objective optimization](./examples/3_basic_multi_objective.cpp)
* [Constrained optimization](./examples/4_constrained_problem.cpp)

-------------------------------------------------------------------------------------------------


