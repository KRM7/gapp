/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_SINGLE_OBJECTIVE_HPP
#define GA_PROBLEMS_SINGLE_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include <vector>
#include <string>
#include <cstddef>

namespace gapp::problems
{
    /**
    * Implementation of the function \f$ -f(\vec{x}) = \langle \vec{x}, \vec{x} \rangle \f$
    * for any number of variables.
    * This is a simple single-objective benchmark function with a single global optimum
    * that is easy to find, and no other local optima.
    * 
    * Evaluated on the hypercube \f$ x_i \in [-5.12, 5.12] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (0, 0, ... , 0) \f$.
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    */
    class Sphere final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a sphere function.
        * 
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Sphere(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Sphere", Bounds{ -5.12, 5.12 }, std::vector(num_vars, 0.0), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %Rastrigin function for any number of variables,
    * modified for maximization.
    * This function has many uniformly distributed local optima and a single global optimum.
    * 
    * \f[ -f(\vec{x}) = 10d + \sum_{i=1}^d \left[ x_i^2 - 10\cos(2\pi x_i) \right] \f]
    * 
    * Evaluated on the hypercube \f$ x_i \in [-5.12, 5.12] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (0, 0, ... , 0) \f$.
    * 
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    * 
    * @see
    *   %Rastrigin, L. A. "Systems of extremal control." Nauka, Moscow (1974).
    */
    class Rastrigin final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %Rastrigin function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Rastrigin(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Rastrigin", Bounds{ -5.12, 5.12 }, std::vector(num_vars, 0.0), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %Rosenbrock function for any number of variables,
    * modified for maximization.
    * This single-objective benchmark function has no local optima, and only a single global optimum,
    * but this optimum can be very difficult to find.
    * 
    * \f[ -f(\vec{x}) = \sum_{i=1}^{d-1} \left[ 100( x_{i+1} - x_i^2 )^2 + (x_i - 1)^2 \right] \f]
    * 
    * Evaluated on the hypercube \f$ x_i \in [-2.048, 2.048] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (0, 0, ... , 0) \f$.
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    * 
    * @see
    *   %Rosenbrock, H. H. "An automatic method for finding the greatest or least value of a function."
    *   The computer journal 3, no. 3 (1960): 175-184.
    */
    class Rosenbrock final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %Rosenbrock function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Rosenbrock(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Rosenbrock", Bounds{ -2.048, 2.048 }, std::vector(num_vars, 1.0), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %Schwefel function for any number of variables,
    * modified for maximization.
    * This single-objective benchmark function has many local optima,
    * and a global optimum that is far away from the next best local optima.
    *
    * \f[ -f(\vec{x}) = 418.98d - \sum_{i=1}^d x_i \sin\left( \sqrt{|x_i|} \right) \f]
    * 
    * Evaluated on the hypercube \f$ x_i \in [-500.0, 500.0] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (420.9687, 420.9687, ... , 420.9687) \f$.
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    */
    class Schwefel final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %Schwefel function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Schwefel(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Schwefel", Bounds{ -500.0, 500.0 }, std::vector(num_vars, 420.9687), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<double>& vars) const override;
    };


    /**
    * Implementation of the %Griewank function for any number of variables,
    * modified for maximization.
    * This single-objective benchmark function has many uniformly distributed local optima
    * and a single global optimum.
    * While it can be used with any number of variables, the difficulty of finding the global
    * optimum does not increase significantly for variable counts above ~10.
    *
    * \f[ -f(\vec{x}) = 1 + \sum_{i=1}^d\frac{x_i^2}{4000} - \prod_{i=1}^d\cos\left( \frac{x_i}{\sqrt{i}} \right) \f]
    * 
    * Evaluated on the hypercube \f$ x_i \in [-600.0, 600.0] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (0, 0, ... , 0) \f$.
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    * 
    * @see
    *   Locatelli, M. "A note on the Griewank test function." Journal of global optimization 25, no. 2 (2003): 169-174.
    * @see
    *   Huang, Y., et al. "Unusual phenomenon of optimizing the Griewank function with the increase of dimension."
    *   Frontiers of Information Technology & Electronic Engineering 20, no. 10 (2019): 1344-1360.
    */
    class Griewank final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %Griewank function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Griewank(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Griewank", Bounds{ -600.0, 600.0 }, std::vector(num_vars, 0.0), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %Ackley function for any number of variables,
    * modified for maximization.
    * This single-objective benchmark function has a single global optimum, and many local optima.
    * The regions further away from the optimum are nearly flat.
    * 
    * \f[ -f(\vec{x}) = 20 + e - 20e^{\left( -0.2\sqrt{ \frac{1}{d}\sum_{i=1}^d x_i^2 } \right)} -
    *                   e^{\left( \frac{1}{d}\sum_{i=1}^d\cos(2\pi x_i) \right)} \f]
    *
    * Evaluated on the hypercube \f$ x_i \in [-32.768, 32.768] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (0, 0, ... , 0) \f$.
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   %Ackley, D H. "A connectionist machine for genetic hillclimbing." (1987).
    */
    class Ackley final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create an %Ackley function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Ackley(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Ackley", Bounds{ -32.768, 32.768 }, std::vector(num_vars, 0.0), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the Lévy function for any number of variables,
    * modified for maximization.
    * 
    * \f[ -f(\vec{x}) = \sin^2(\pi w_1) + (w_d - 1)^2( 1 + \sin^2(2\pi w_d) ) + \sum_{i=1}^d (w_i - 1)^2 (1 + 10\sin^2(\pi w_i + 1)) \f]
    * \f[ \textrm{where  } w_i = 1 + \frac{x_i - 1}{4} \f]
    *
    * Evaluated on the hypercube \f$ x_i \in [-10.0, 10.0] \f$.
    * 
    * The global optimum is \f$ f(\vec{x}) = 0 \f$ at \f$ \vec{x} = (1, 1, ... , 1) \f$.
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    */
    class Levy final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a Lévy function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Levy(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunction("Levy", Bounds{ -10.0, 10.0 }, std::vector(num_vars, 1.0), 0.0, bits_per_var)
        {}

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace gapp::problems

#endif // !GA_PROBLEMS_SINGLE_OBJECTIVE_HPP