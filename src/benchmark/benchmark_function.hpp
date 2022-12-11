/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BENCHMARK_BENCHMARK_FUNCTION_HPP
#define GA_BENCHMARK_BENCHMARK_FUNCTION_HPP

#include "../core/ga_base.decl.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/utility.hpp"
#include <vector>
#include <string>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::benchmark
{
    template<typename T>
    class BenchmarkFunction
    {
    public:
        using Bounds = GA<T>::GeneBounds;

        BenchmarkFunction(std::string name, size_t num_objs, size_t num_vars, const std::vector<Bounds>& bounds) :
            name_(std::move(name)), bounds_(bounds), num_objs_(num_objs), num_vars_(num_vars)
        {
            if (num_objs == 0) GA_THROW(std::invalid_argument, "Number of objectives must be at least 1.");
            if (num_vars == 0) GA_THROW(std::invalid_argument, "Number of variables must be at least 1.");
            if (bounds.size() != num_vars) GA_THROW(std::invalid_argument, "Mismatching number of variables and bounds vector sizes.");
        }

        BenchmarkFunction(std::string name, size_t num_objs, size_t num_vars, Bounds bounds) :
            BenchmarkFunction(std::move(name), num_objs, num_vars, std::vector(num_vars, bounds))
        {}

        virtual ~BenchmarkFunction() = default;

        const std::vector<Bounds>& bounds() const noexcept { return bounds_; }

        size_t num_obj() const noexcept  { return num_objs_; }
        size_t num_vars() const noexcept { return num_vars_; }

        const std::string& name() const { return name_; }

        std::vector<double> operator()(const std::vector<T>& x) const
        {
            assert(x.size() == num_vars_);
            return invoke(x);
        }

    protected:
        virtual std::vector<double> invoke(const std::vector<T>& x) const = 0;

        std::vector<Bounds> bounds_;

    private:
        std::string name_;
        size_t num_objs_;
        size_t num_vars_;
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_BENCHMARK_FUNCTION_HPP