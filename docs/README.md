# Genetic Algorithms in C++ 

<p> Implementations of genetic algorithms in C++. </p>

## Overview

<p>
This library implements several types of genetic algorithms in C++ for solving both single-
and multi-objective optimization problems. Genetic algorithms with several commonly used
encoding methods are implemented, along with commonly used selection, crossover, and mutation
methods for each encoding method. The user can either choose one of these already implemented 
methods to use, or they can define their own genetic operators.
</p>

<p> 
Most parts of the algorithms can be customized to use user-defined methods/types,
including the definition of custom, user-defined encoding types,
custom, user-defined selection methods for the single-objective algorithms,
and user-defined crossover and mutation methods.
A repair function can also be defined in order to create memetic algorithms.
</p>

<p>
The genetic algorithms are parallelized through the use of C++17 parallel algorithms.
</p>

#### Requirements

<p>
The library consists of header files only, and does not depend on any other libraries.
The only requirement is a compiler that supports the C++20 standard.
</p>

## Usage examples

The [examples](/examples) folder contains some simple usage examples for 
every major feature implemented in the library.

## Encoding types

<p>
The base GA class can be inherited from to create genetic algorithms with any encoding type,
but the encoding types already implemented are:
</p>

* Binary encoding
* Real encoding
* Permutation encoding
* Integer encoding

<p>
Each of these differently encoded genetic algorithms are implemented as a separate class
fully implemented in their respective header files. <br/>
The genetic algorithms are used by instantiating one of these classes.
<p/>


## Available algorithms
Both single- and multi-objective algorithms are available and can be used regardless
of the encoding. These algorithms are:
* Single-objective elitist genetic algorithm
* Non-Dominated Sorting Genetic Algorithm II (NSGA-II)
* Non-Dominated Sorting Genetic Algorithm III (NSGA-III)


## Genetic operators

### Crossover and mutation methods
Custom crossover and mutation functions implemented by the user can be used, but
for each encoding type, some commonly used crossover and mutation methods are already implemented:

|                    | Crossover methods | Mutation methods |
|:------------------:|:-----------------:|:----------------:|
| **Binary GA**      |*-single-point*<br/>*-two-point*<br/>*-n-point*<br/>*-uniform*    |*-standard*|
| **Real GA**        |*-arithmetic*<br/>*-BLX-α*<br/>*-simulated binary*<br/>*-Wright*  |*-random*<br/>*-polynomial*<br/>*-gauss*<br/>*-non-uniform*<br/>*-boundary*|
| **Permutation GA** |*-order (OX1)*<br/>*-cycle (CX)*<br/>*-edge assembly (EAX)*<br/>*-partially mapped (PMX)* |*-single-swap*<br/>*-scramble*<br/>*-inversion*|
| **Integer GA**     |*-single-point*<br/>*-two-point*<br/>*-n-point*<br/>*-uniform*    |*-standard*|


### Selection methods (SOGA)
<p>
The multi-objective algorithms (NSGA-II and NSGA-III) always use their pre-defined
selection methods and they can not be changed, 
but the selection method in the single-objective algorithm can either be a user implemented method
or one of the following already implemented methods:
</p>

 * Roulette-wheel selection
 * Tournament selection
 * Linear rank selection
 * Sigma-scaling selection
 * Boltzmann selection

### Other

The stop condition used in the algorithms can be chosen by the user from a set of implemented
stop conditions, and the initial population can also be set by the user instead of being
randomly generated. See the [examples](/examples).


## References
<p>NSGA-II:</p>

* Deb, K., et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II." 
*IEEE transactions on evolutionary computation* 6.2 (2002): 182-197.

<p>NSGA-III:</p>

* Deb, K., and Jain, H. "An evolutionary many-objective optimization algorithm using reference-point-based nondominated sorting approach, part I: solving problems with box constraints." 
*IEEE transactions on evolutionary computation* 18.4 (2013): 577-601.

* Blank, J., Deb, K., & Roy, P. C. "Investigating the normalization procedure of NSGA-III."
*International Conference on Evolutionary Multi-Criterion Optimization* (2019): 229-240.

* Seada, H., & Deb, K. "U-NSGA-III: a unified evolutionary optimization procedure for single, multiple, and many objectives: proof-of-principle results."
*International conference on evolutionary multi-criterion optimization* (2015): 34-49.

<p>Other:</p>

* Deb, K., Bandaru, S., & Seada, H. "Generating uniformly distributed points on a unit simplex for evolutionary many-objective optimization."
*In International Conference on Evolutionary Multi-Criterion Optimization* (2019): 179-190.

* Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
*Evolutionary multiobjective optimization* (2005): 105-145.

* He, Linjun, et al. "Dynamic Normalization in MOEA/D for Multiobjective optimization." 
*IEEE Congress on Evolutionary Computation (CEC)* (2020): 1-8.
