/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_MULTI_OBJECTIVE_HPP
#define GA_PROBLEMS_MULTI_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::problems
{
    /**
    * Implementation of the %Kursawe function for any number of variables,
    * modified for maximization. It has 2 objectives, and the pareto front
    * is made up of multiple disconnected segments.
    *
    * Evaluated on the hypercube \f$ x_i \in [-5.0,\ 5.0] \f$.
    *
    * The approximate extreme points in the objective-space are:
    *    \f[ \textrm{ ideal-point: } (  10(n_{vars} - 1),\ 3.85(n_{vars} - 1) + 4)  \f]
    *    \f[ \textrm{ nadir-point: } (7.25(n_{vars} - 1),\           0.0         )  \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    * 
    * @see
    *   %Kursawe, F. "A variant of evolution strategies for vector optimization."
    *   International conference on parallel problem solving from nature (1991): 193-197
    */
    class Kursawe final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %Kursawe function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit Kursawe(size_t num_vars = 3, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %ZDT1 function for any number of variables,
    * modified for maximization. This 2-objective benchmark function has a continuous,
    * convex pareto front in the objective-space.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0,\ 1.0] \f$.
    * 
    * The optimal solutions are \f$ x_1 \in [0.0,\ 1.0] \f$, and \f$ x_{rest} = 0.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT1 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %ZDT1 function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit ZDT1(size_t num_vars = 30, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %ZDT2 function for any number of variables,
    * modified for maximization. This 2-objective benchmark function has a continuous,
    * non-convex pareto front in the objective-space.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0,\ 1.0] \f$.
    * 
    * The optimal solutions are \f$ x_1 \in [0.0,\ 1.0] \f$, and \f$ x_{rest} = 0.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT2 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %ZDT2 function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit ZDT2(size_t num_vars = 30, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %ZDT3 function for any number of variables,
    * modified for maximization. This 2-objective benchmark function has a
    * discontinuous pareto front made up of 5 disconnected segments in the objective-space.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0,\ 1.0] \f$.
    * 
    * The optimal solutions are \f$ x_1 \in [0.0,\ 0.083] \cup [0.182,\ 0.258] \cup
    * [0.409,\ 0.454] \cup [0.618\, 0.653] \cup [0.823,\ 0.852] \f$, and \f$ x_{rest} = 0.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.00,\ 0.8) \f]
    *   \f[ \textrm{ nadir-point: } (-0.85, -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT3 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %ZDT3 function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit ZDT3(size_t num_vars = 30, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %ZDT4 function for any number of variables,
    * modified for maximization. This 2-objective benchmark function is
    * the most difficult problem of the ZDT suite. It has a large number of local
    * pareto fronts in the objective space, which makes it easy for the %GAs to get
    * stuck along one of these local fronts.
    *
    * Evaluated on \f$ x_1 \in [0.0,\ 1.0] \f$, and \f$ x_{rest} \in [-5.0,\ 5.0] \f$.
    * 
    * The optimal solutions are \f$ x_1 \in [0.0,\ 1.0] \f$, and \f$ x_{rest} = 0.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT4 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %ZDT4 function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit ZDT4(size_t num_vars = 10, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };


    /**
    * Implementation of the %ZDT5 function for any number of variables,
    * modified for maximization. In this 2-objective benchmark function,
    * the variables are binary strings, not real numbers like they are in the rest
    * of the functions of the ZDT suite.
    *
    * The optimal solutions are \f$ x_1 = anything \f$, and \f$ x_{rest} = ones \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( -1.0, -\frac{(n_{vars} - 1)}{31})  \f]
    *   \f[ \textrm{ nadir-point: } (-31.0,        -n_{vars} + 1      )  \f]
    *
    * This benchmark function can only be used with the binary-encoded %GA,
    * unlike the rest of the benchmark functions in the ZDT suite.
    *
    * @see
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT5 final : public BenchmarkFunction<BinaryGene>
    {
    public:
        /**
        * Create a %ZDT5 function.
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
    * Implementation of the %ZDT6 function for any number of variables,
    * modified for maximization. This 2-objective benchmark function has a
    * non-convex pareto front in the objective space, along which the solutions
    * are distributed non-uniformly.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are \f$ x_1 \in [0.0, 1.0] \f$, and \f$ x_{rest} = 0.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0 ) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -0.92) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Zitzler, E., Deb, K., and Thiele, L. "Comparison of multiobjective evolutionary algorithms: Empirical results."
    *   Evolutionary computation 8, no. 2 (2000): 173-195.
    */
    class ZDT6 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Create a %ZDT6 function.
        *
        * @param num_vars The number of variables. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit ZDT6(size_t num_vars = 10, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_MULTI_OBJECTIVE_HPP