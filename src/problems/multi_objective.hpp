/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/**
* Implementations of some commonly used multi-objective benchmark functions
* that can be used for testing the binary- and real-encoded algorithms. \n
*
* All of the functions are implemented for maximization, and can be constructed for any
* number of variables, but they all have 2 objectives. \n
* Most of them can be used for both the binary- and real-encoded algorithms.
*/

#ifndef GA_PROBLEMS_MULTI_OBJECTIVE_HPP
#define GA_PROBLEMS_MULTI_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::problems
{
    /**
    * Implementation of the Kursawe function for any number of variables (and 2 objectives). \n
    *
    * Evaluated on the hypercube x_i = [-5.0, 5.0]. \n
    * The function is implemented for maximization. It has multiple disconnected pareto fronts. \n
    * The approximate extreme points in the objective-space are: \n
    *   ideal-point: [10*(nvars-1), 3.85(nvars-1) + 4] \n
    *   nadir-point: [7.25*(nvars-1), 0.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    * 
    * See:
    *   Kursawe, F. "A variant of evolution strategies for vector optimization."
    *   International conference on parallel problem solving from nature (1991): 193-197
    */
    class Kursawe final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a Kursawe function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit Kursawe(size_t num_vars = 3, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT1 function for any number of variables (and 2 objectives). \n
    * The function has a contiguous, convex pareto front in the objective-space. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT1 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a ZDT1 objective function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT1(size_t num_vars = 30, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT2 function for any number of variables (and 2 objectives). \n
    * The function has a contiguous, non-convex pareto front in the objective-space. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT2 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a ZDT2 objective function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT2(size_t num_vars = 30, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT3 function for any number of variables (and 2 objectives). \n
    * The function has a discontinuous pareto front of 5 fronts in the objective-space. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are
    * x_1 = [0.0, 0.0830] U [0.1822, 0.2577] U [0.4093, 0.4538] U [0.6183, 0.6525] U [0.8233, 0.8518],
    * and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.00,  0.8] \n
    *   nadir-point: [-0.85, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT3 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a ZDT3 objective function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT3(size_t num_vars = 30, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT4 function for any number of variables (and 2 objectives). \n
    * This is the most difficult problem of the ZDT suite,
    * the function has a large number of local pareto fronts in the objective space,
    * which makes it easy for the algorithms to get stuck along one of these local fronts. \n
    *
    * Evaluated on x_1 = [0.0, 1.0] and x_rest = [-5.0, 5.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -1.0] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT4 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a ZDT4 objective function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT4(size_t num_vars = 10, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the ZDT5 function for any number of variables (and 2 objectives). \n
    * Unlike the rest of the functions in the ZDT suite, the variables are binary strings, not
    * real numbers. \n
    *
    * The function is implemented for maximization, the optimal solutions are x_1 = anything, and x_rest = all_ones (binary chromosomes). \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ -1.0, -(nvars-1)/31] \n
    *   nadir-point: [-31.0,    -nvars+1  ] \n
    *
    * This benchmark function can only be used with the binary-encoded multi-objective algorithms,
    * unlike the rest of the benchmark functions in the ZDT suite. \n
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT5 final : public BenchmarkFunction<BinaryGene>
    {
    public:
        /**
        * Construct a ZDT5 objective function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        */
        explicit ZDT5(size_t num_vars = 11);

    private:
        FitnessVector invoke(const std::vector<BinaryGene>& vars) const override;

        static constexpr size_t FIRST_BITS = 30;
        static constexpr size_t REST_BITS = 5;
    };


    /**
    * Implementation of the ZDT6 function for any number of variables (and 2 objectives). \n
    * The function has a nonconvex pareto front in the objective space,
    * along which the solutions are non-uniformly distributed. \n
    *
    * Evaluated on x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are x_1 = [0.0, 1.0] and x_rest = 0.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: [ 0.0,  0.0] \n
    *   nadir-point: [-1.0, -0.92] \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT6 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a ZDT6 objective function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit ZDT6(size_t num_vars = 10, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_MULTI_OBJECTIVE_HPP