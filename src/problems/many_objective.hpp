/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_MANY_OBJECTIVE_HPP
#define GA_PROBLEMS_MANY_OBJECTIVE_HPP

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include <cstddef>

namespace gapp::problems
{
    /**
    * Implementation of the %DTLZ1 function for any number of objectives, modified
    * for maximization. This is the simplest test problem of the DTLZ test suite,
    * with a linear pareto front and a large number (\f$ 11^K \f$) of local pareto fronts.
    * 
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.5 \textrm{ for } i >= n_{obj}-1 \f]
    * The pareto optimal front in the objective space is the linear hyperplane
    * where \f$ \sum_m f_m = 0.5 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ...,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-0.5, -0.5,\ ..., -0.5) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ1 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ1 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ1(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 5;
    };


    /**
    * Implementation of the %DTLZ2 function for any number of objectives, modified
    * for maximization. The problem has a non-linear, spherical pareto front.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are: \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.5 \textrm{ for } i >= n_{obj}-1 \f]
    * The optimal front in the objective space is the surface of the unit hypersphere:
    * \f$ \sum_m f_m^2 = 1.0. \f$
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ...,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0,\ ..., -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ2 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ2 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ2(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the %DTLZ3 function for any number of objectives, modified
    * for maximization. This problem is a modified version of the DTLZ2 problem with
    * many local pareto-optimal fronts, making it more difficult to find the global optimal front.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are: \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.5 \textrm{ for } i >= n_{obj}-1 \f]
    * The pareto optimal front in the objective space is the surface of the unit hypersphere:
    * \f$ \sum_m f_m^2 = 1.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ...,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0,\ ..., -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded multi-objective algorithms. \n
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ3 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ3 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ3(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the %DTLZ4 function for any number of objectives, modified
    * for maximization. This problem is a modification of the DTLZ2 problem,
    * where the solutions are not uniformly distributed along the pareto-optimal front.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are: \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.5 \textrm{ for } i >= n_{obj}-1 \f]
    * The optimal front in the objective space is the surface of the unit hypersphere:
    * \f$ \sum_m f_m^2 = 1.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ...,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0,\ ..., -1.0) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ4 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ4 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ4(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the %DTLZ5 function for any number of objectives, modified
    * for maximization. This problem is a modified version of the DTLZ2 problem,
    * where the pareto-optimal front is a degenerated curve instead of the surface of a hypersphere.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are: \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.5 \textrm{ for } i >= n_{obj}-1 \f]
    * The pareto optimal solutions are along a curve that satisfies \f$ \sum_m f_m^2 = 1.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ...,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } \left( \left( \frac{-1.0}{\sqrt{2}} \right)^{n_{obj}-2},
    *                                      \left( \frac{-1.0}{\sqrt{2}} \right)^{n_{obj}-2},
    *                                      \left( \frac{-1.0}{\sqrt{2}} \right)^{n_{obj}-3},
    *                                      \ ...,
    *                                       -1.0 \right) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ5 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ5 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ5(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the %DTLZ6 function for any number of objectives, modified
    * for maximization. This problem is a modified version of the DTLZ5 problem with
    * many local pareto-optimal fronts, making it more difficult to find the global one.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are: \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.0 \textrm{ for } i >= n_{obj}-1 \f]
    * The pareto optimal solutions are along a curve that satisfies \f$ \sum_m f_m^2 = 1.0 \f$.
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ...,\ 0.0) \f]
    *   \f[ \textrm{ nadir-point: } \left( \left( \frac{-1.0}{\sqrt{2}} \right)^{n_{obj}-2},
    *                                      \left( \frac{-1.0}{\sqrt{2}} \right)^{n_{obj}-2},
    *                                      \left( \frac{-1.0}{\sqrt{2}} \right)^{n_{obj}-3},
    *                                      \ ...,
    *                                       -1.0 \right) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ6 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ6 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ6(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 10;
    };


    /**
    * Implementation of the %DTLZ7 function for any number of objectives, modified
    * for maximization. This problem has (\f$ 2^{n_{obj}-1} \f$) disconnected
    * pareto fronts in the objective space instead of just a single continuous one.
    *
    * Evaluated on the hypercube \f$ x_i \in [0.0, 1.0] \f$.
    * 
    * The optimal solutions are: \f[ x_i \in [0.0, 1.0] \textrm{ for } 0 <= i < n_{obj}-1,
    * \textrm{ and } x_i = 0.0 \textrm{ for } i >= n_{obj}-1 \f]
    * 
    * The extreme points in the objective-space are:
    *   \f[ \textrm{ ideal-point: } ( 0.0,\ 0.0,\ ..., -0.307n_{obj}-1.693) \f]
    *   \f[ \textrm{ nadir-point: } (-1.0, -1.0,\ ...,       -2n_{obj}    ) \f]
    *
    * This benchmark function can be used for both the real- and binary-encoded %GAs.
    *
    * @see
    *   Deb, K., et al. "Scalable test problems for evolutionary multiobjective optimization."
    *   Evolutionary multiobjective optimization (2005), pp. 105-145.
    * @see
    *   Deb, K., et al. "Scalable multi-objective optimization test problems."
    *   Proceedings of the 2002 Congress on Evolutionary Computation. vol. 1, pp. 825-830.
    */
    class DTLZ7 final : public BenchmarkFunction<RealGene>
    {
    public:
        /**
        * Construct a %DTLZ7 objective function.
        *
        * @param num_obj The number of objectives. Must be at least 2.
        * @param bits_per_var The number of bits representing a variable when used with the binary-encoded %GA.
        */
        explicit DTLZ7(size_t num_obj, size_t bits_per_var = 32);

    private:
        FitnessVector invoke(const Chromosome<RealGene>& vars) const override;

        static constexpr size_t K = 20;
    };

} // namespace gapp::problems

#endif // !GA_PROBLEMS_MANY_OBJECTIVE_HPP