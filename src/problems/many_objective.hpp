/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/**
* Implementations of some commonly used multi/many-objective benchmark functions
* that can be used for testing the binary- and real-encoded algorithms. \n
*
* All of the functions are implemented for maximization, and can be constructed for any
* number of objectives. The number of variables used is determined by the nuumber of objectives. \n
* They can be used for both the binary- and real-encoded algorithms.
*/

#ifndef GA_PROBLEMS_MANY_OBJECTIVE_HPP
#define GA_PROBLEMS_MANY_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include <vector>
#include <cstddef>

namespace genetic_algorithm::problems
{
    /**
    * Implementation of the DTLZ1 function for any number of objectives. \n
    * This is the simplest test problem of the DTLZ test suite, with a linear pareto front
    * and a large number (11^K) of local pareto-optimal fronts. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.5 for i >= nobj-1,
    * and the pareto optimal front in the objective space is the linear hyperplane where sum(f_m) = 0.5. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ...,  0.0) \n
    *   nadir-point: (-0.5, -0.5, ..., -0.5) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * 
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ1 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ1 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ1(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 5;
    };


    /**
    * Implementation of the DTLZ2 function for any number of objectives. \n
    * The problem has a non-linear, spherical pareto-optimal front. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.5 for i >= nobj-1,
    * and the pareto optimal front in the objective space is the surface of the unit hypersphere: sum(f_m^2) = 1.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ...,  0.0) \n
    *   nadir-point: (-1.0, -1.0, ..., -1.0) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    *
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ2 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ2 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ2(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the DTLZ3 function for any number of objectives. \n
    * This problem is a modified version of the DTLZ2 problem with many local
    * pareto-optimal fronts, making it more difficult to find the global optimal front. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.5 for i >= nobj-1,
    * and the pareto optimal front in the objective space is the surface of the unit hypersphere: sum(f_m^2) = 1.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ...,  0.0) \n
    *   nadir-point: (-1.0, -1.0, ..., -1.0) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    *
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ3 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ3 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ3(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the DTLZ4 function for any number of objectives. \n
    * This problem is a modification of the DTLZ2 problem, where the solutions
    * are not uniformly distributed along the pareto-optimal front. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.5 for i >= nobj-1,
    * and the pareto optimal front in the objective space is the surface of the unit hypersphere: sum(f_m^2) = 1.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ...,  0.0) \n
    *   nadir-point: (-1.0, -1.0, ..., -1.0) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    *
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ4 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ4 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ4(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the DTLZ5 function for any number of objectives. \n
    * This problem is a modified version of the DTLZ2 problem, where the pareto-optimal
    * front is a degenerated curve instead of the surface of a hypersphere. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.5 for i >= nobj-1,
    * and the pareto optimal solutions satisfy: sum(f_m^2) = 1.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ...,  0.0) \n
    *   nadir-point: ( (-1.0/sqrt(2))^(nobj-2), (-1.0/sqrt(2))^(nobj-2), (-1.0/sqrt(2))^(nobj-3), ..., -1.0) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    *
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ5 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ5 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ5(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the DTLZ6 function for any number of objectives. \n
    * This problem is a modified version of the DTLZ5 problem with many local
    * pareto-optimal fronts, making it more difficult to find the global one.
    * (The optimal front is still just a degenerated curve.) \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.0 for i >= nobj-1,
    * and the pareto optimal solutions satisfy: sum(f_m^2) = 1.0. \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ...,  0.0) \n
    *   nadir-point: ( (-1.0/sqrt(2))^(nobj-2), (-1.0/sqrt(2))^(nobj-2), (-1.0/sqrt(2))^(nobj-3), ..., -1.0) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    *
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ6 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ6 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ6(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the DTLZ7 function for any number of objectives. \n
    * This problem has multiple ( 2^(nobj-1) ) disconnected pareto-optimal fronts instead of a single
    * contiguous one. \n
    *
    * Evaluated on the hypercube x_i = [0.0, 1.0]. \n
    * The function is implemented for maximization, the optimal solutions are:
    * x_i = [0.0, 1.0] for 0 <= i < nobj-1, and x_i = 0.0 for i >= nobj-1 \n
    * The extreme points in the objective-space are: \n
    *   ideal-point: ( 0.0,  0.0, ..., -0.307*nobj-1.693) \n
    *   nadir-point: (-1.0, -1.0, ..., -2*nobj) \n
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * See:
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    *
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ7 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a DTLZ7 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded algorithms.
        */
        explicit DTLZ7(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const std::vector<RealGene>& vars) const override;

        static constexpr size_t K = 20;
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_MANY_OBJECTIVE_HPP