/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_PROBLEMS_BENCHMARK_FUNCTION_HPP
#define GA_PROBLEMS_BENCHMARK_FUNCTION_HPP

#include "../core/ga_base.decl.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/math.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <string>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
{

    template<typename T>
    class BenchmarkFunction
    {
    public:
        using Gene      = T;
        using Bounds    = typename GA<T>::GeneBounds;
        using BoundsVec = std::vector<Bounds>;
        using Point     = math::Point;

        BenchmarkFunction(std::string name, size_t num_objs, const BoundsVec& bounds) :
            bounds_(bounds), name_(std::move(name)), num_objs_(num_objs), num_vars_(bounds.size())
        {
            if (num_objs_ == 0) GA_THROW(std::invalid_argument, "Number of objectives must be at least 1.");
            if (num_vars_ == 0) GA_THROW(std::invalid_argument, "Number of variables must be at least 1.");
        }

        BenchmarkFunction(std::string name, size_t num_objs, size_t num_vars, Bounds bounds) :
            BenchmarkFunction(std::move(name), num_objs, std::vector(num_vars, bounds))
        {}

        virtual ~BenchmarkFunction() = default;

        const BoundsVec& bounds() const noexcept { return bounds_; }

        size_t num_obj() const noexcept  { return num_objs_; }
        size_t num_vars() const noexcept { return num_vars_; }

        const std::string& name() const { return name_; }

        const Point& optimal_value() const noexcept { return optimal_value_; }
        const std::vector<Gene>& optimum() const noexcept { return optimum_; }

        std::vector<double> operator()(const std::vector<T>& x) const
        {
            assert(x.size() == num_vars_);
            return invoke(x);
        }

    protected:
        virtual std::vector<double> invoke(const std::vector<T>& x) const = 0;

        std::vector<Bounds> bounds_;
        std::vector<Gene> optimum_;
        Point optimal_value_;

    private:
        std::string name_;
        size_t num_objs_;
        size_t num_vars_;
    };


    /** ... */
    class BenchmarkFunctionReal1 : public BenchmarkFunction<RealGene>
    {
    public:
        BenchmarkFunctionReal1(std::string name, const std::vector<Bounds>& bounds, size_t bits_per_var, double optimal_value, const std::vector<Gene>& optimum) :
            BenchmarkFunction<Gene>(std::move(name), 1, bounds), var_bits_(bits_per_var)
        {
            if (optimum.size() != num_vars()) GA_THROW(std::invalid_argument, "Mismatching number of variables and optimum vector sizes.");

            optimum_ = optimum;
            optimal_value_ = Point(1, optimal_value);
        }

        BenchmarkFunctionReal1(std::string name, size_t num_vars, Bounds bounds, size_t bits_per_var, double optimal_value, const std::vector<Gene>& optimum) :
            BenchmarkFunctionReal1(std::move(name), std::vector(num_vars, bounds), bits_per_var, optimal_value, optimum)
        {}

        size_t num_bits() const noexcept { return num_vars() * var_bits_; }
        size_t var_bits() const noexcept { return var_bits_; }

        using BenchmarkFunction<Gene>::operator();
        std::vector<double> operator()(const std::vector<BinaryGene>& binary_chrom) const { return invoke(convert(binary_chrom)); }

    private:
        std::vector<RealGene> convert(const std::vector<BinaryGene>& binary_chrom) const;

        size_t var_bits_;
    };


    /** ... */
    class BenchmarkFunctionRealN : public BenchmarkFunction<RealGene>
    {
    public:
        BenchmarkFunctionRealN(std::string name, size_t num_obj, const std::vector<Bounds>& bounds, size_t bits_per_var) :
            BenchmarkFunction<RealGene>(std::move(name), num_obj, bounds), var_bits_(bits_per_var)
        {
            if (num_obj < 2) GA_THROW(std::invalid_argument, "Not enough objectives for a multi-objective benchmark functions.");
        }

        BenchmarkFunctionRealN(std::string name, size_t num_obj, size_t num_vars, Bounds bounds, size_t bits_per_var) :
            BenchmarkFunctionRealN(std::move(name), num_obj, std::vector(num_vars, bounds), bits_per_var)
        {}

        size_t num_bits() const noexcept { return num_vars() * var_bits_; }
        size_t var_bits() const noexcept { return var_bits_; }

        const Point& ideal_point() const noexcept { return ideal_point_; }
        const Point& nadir_point() const noexcept { return nadir_point_; }

        using BenchmarkFunction<Gene>::operator();
        std::vector<double> operator()(const std::vector<BinaryGene>& binary_chrom) const { return invoke(convert(binary_chrom)); }

    protected:
        std::vector<RealGene> convert(const std::vector<BinaryGene>& binary_chrom) const;

        Point ideal_point_;
        Point nadir_point_;
        size_t var_bits_;
    };


    /** ... */
    class BenchmarkFunctionBinaryN : public BenchmarkFunction<BinaryGene>
    {
    public:
        BenchmarkFunctionBinaryN(std::string name, size_t num_obj, size_t num_bits) :
            BenchmarkFunction<BinaryGene>(std::move(name), num_obj, num_bits, Bounds{ 0, 1 })
        {}

        const Point& ideal_point() const noexcept { return ideal_point_; }
        const Point& nadir_point() const noexcept { return nadir_point_; }

    protected:
        Point ideal_point_;
        Point nadir_point_;
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_BENCHMARK_FUNCTION_HPP