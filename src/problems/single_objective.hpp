/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/**
* Implementations of some commonly used single-objective benchmark functions
* that can be used for testing the binary- and real-encoded algorithms. \n
* 
* All of the functions are implemented for maximization, and can be constructed for any
* number of variables. \n
* They can be used for both the binary- and real-encoded algorithms.
*/

#ifndef GA_PROBLEMS_SINGLE_OBJECTIVE_HPP
#define GA_PROBLEMS_SINGLE_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <utility>
#include <string>
#include <cstddef>

namespace genetic_algorithm::problems
{
    /** ... */
    class BenchmarkFunctionReal1 : public BenchmarkFunction<RealGene>
    {
    public:
        explicit BenchmarkFunctionReal1(std::string name, size_t num_vars, Bounds bounds, size_t bits_per_var, double optimal_value, const std::vector<RealGene>& optimum) :
            BenchmarkFunction<RealGene>(std::move(name), 1, num_vars, bounds), optimum_(optimum), optimal_value_(1, optimal_value), var_bits_(bits_per_var)
        {
            if (optimum.size() != num_vars) GA_THROW(std::invalid_argument, "Mismatching number of variables and optimum vector sizes.");
        }

        explicit BenchmarkFunctionReal1(std::string name, size_t num_vars, const std::vector<Bounds>& bounds, size_t bits_per_var, double optimal_value, const std::vector<RealGene>& optimum) :
            BenchmarkFunction<RealGene>(std::move(name), 1, num_vars, bounds), optimum_(optimum), optimal_value_(1, optimal_value), var_bits_(bits_per_var)
        {
            if (optimum.size() != num_vars) GA_THROW(std::invalid_argument, "Mismatching number of variables and optimum vector sizes.");
        }

        const std::vector<RealGene>& optimal_value() const noexcept { return optimal_value_; }
        const std::vector<RealGene>& optimum() const noexcept { return optimum_; }

        size_t num_bits() const noexcept { return num_vars() * var_bits_; }
        size_t var_bits() const noexcept { return var_bits_; }

        using BenchmarkFunction<double>::operator();
        std::vector<double> operator()(const std::vector<BinaryGene>& binary_chrom) const { return invoke(convert(binary_chrom)); }

    private:
        std::vector<double> convert(const std::vector<BinaryGene>& binary_chrom) const;

        std::vector<RealGene> optimum_;
        std::vector<RealGene> optimal_value_;
        size_t var_bits_;
    };


    /**
    * Implementation of a simple sphere function (x^2) for any number of variables. \n
    * This is a simple benchmark function with a single global optimum that is easy to find, and
    * no other local optima. \n
    * 
    * Evaluated on the hypercube x_i = [-5.12, 5.12]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (0, 0, ... , 0). \n
    *
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    */
    class Sphere : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct a sphere function.
        * 
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Sphere(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Sphere", num_vars, Bounds{ -5.12, 5.12 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the Rastrigin function for any number of variables. \n
    * The function has many uniformly distributed local optima and a single global optimum. \n
    * 
    * Evaluated on the hypercube x_i = [-5.12, 5.12]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (0, 0, ... , 0). \n
    * 
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    * 
    * See:
    *   Rastrigin, L. A. "Systems of extremal control." Nauka, Moscow (1974).
    */
    class Rastrigin : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct a Rastrigin function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Rastrigin(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Rastrigin", num_vars, Bounds{ -5.12, 5.12 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the Rosenbrock function for any number of variables. \n
    * The function has no local optima, and only a single global optimum, but this
    * optimum can be very difficult to find. \n
    * 
    * Evaluated on the hypercube x_i = [-2.048, 2.048]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (0, 0, ... , 0). \n
    *
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    * 
    * See:
    *   Rosenbrock, H. H. "An automatic method for finding the greatest or least value of a function."
    *   The computer journal 3, no. 3 (1960): 175-184.
    */
    class Rosenbrock : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct a Rosenbrock function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Rosenbrock(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Rosenbrock", num_vars, Bounds{ -2.048, 2.048 }, bits_per_var, 0.0, std::vector(num_vars, 1.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the Schwefel function for any number of variables. \n
    * The function has many local optima, and a global optimum that is far away from the next best local optima. \n
    *
    * Evaluated on the hypercube x_i = [-500.0, 500.0]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (420.9687, 420.9687, ... , 420.9687). \n
    *
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    */
    class Schwefel : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct a Schwefel function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Schwefel(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Schwefel", num_vars, Bounds{ -500.0, 500.0 }, bits_per_var, 0.0, std::vector(num_vars, 420.9687))
        {}

    private:
        std::vector<double> invoke(const std::vector<double>& vars) const override;
    };


    /**
    * Implementation of the Griewank function for any number of variables. \n
    * The function has many uniformly distributed local optima and a single global optimum.
    * While it can be used with any number of variables, the difficulty of finding the global
    * optimum does not increase significantly for variable counts above ~10. \n
    *
    * Evaluated on the hypercube x_i = [-600.0, 600.0]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (0, 0, ... , 0). \n
    *
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    * 
    * See:
    *   Locatelli, M. "A note on the Griewank test function." Journal of global optimization 25, no. 2 (2003): 169-174.
    * 
    *   Huang, Y., et al. "Unusual phenomenon of optimizing the Griewank function with the increase of dimension."
    *   Frontiers of Information Technology & Electronic Engineering 20, no. 10 (2019): 1344-1360.
    */
    class Griewank : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct a Griewank function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Griewank(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Griewank", num_vars, Bounds{ -600.0, 600.0 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the Ackley function for any number of variables. \n
    * The function has a single global optimum, and many local optima. The regions further away
    * from the optimum are nearly flat. \n
    *
    * Evaluated on the hypercube x_i = [-32.768, 32.768]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (0, 0, ... , 0). \n
    *
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    *
    * See:
    *   Ackley, D H. "A connectionist machine for genetic hillclimbing." (1987).
    */
    class Ackley : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct an Ackley function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Ackley(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Ackley", num_vars, Bounds{ -32.768, 32.768 }, bits_per_var, 0.0, std::vector(num_vars, 0.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the Lévy function for any number of variables. \n
    *
    * Evaluated on the hypercube x_i = [-10.0, 10.0]. \n
    * The function is implemented for maximization, the global optimum is f(x) = 0, at x = (1, 1, ... , 1). \n
    *
    * This benchmark function can be used for both the real- and binary-encoded single-objective algorithms.
    */
    class Levy : public BenchmarkFunctionReal1
    {
    public:
        /**
        * Construct a Lévy function.
        *
        * @param num_vars The number of variables. Must be at least 1.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Levy(size_t num_vars, size_t bits_per_var = 32) :
            BenchmarkFunctionReal1("Levy", num_vars, Bounds{ -10.0, 10.0 }, bits_per_var, 0.0, std::vector(num_vars, 1.0))
        {}

    private:
        std::vector<double> invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_SINGLE_OBJECTIVE_HPP