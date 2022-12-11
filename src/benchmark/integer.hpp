/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_BENCHMARK_INTEGER_HPP
#define GA_BENCHMARK_INTEGER_HPP

#include "benchmark_function.hpp"
#include <vector>
#include <string>
#include <utility>
#include <cstddef>
#include <cassert>

/* Benchmark functions for the IntegerGA. */
namespace genetic_algorithm::benchmark
{
    class StringFinder : public BenchmarkFunction<IntegerGene>
    {
    public:
        explicit StringFinder(std::string target) :
            BenchmarkFunction<IntegerGene>("StringFinder", 1, target.size(), Bounds{ 0, 95 }), target_(std::move(target))
        {}

        double optimal_value() const noexcept { return double(num_vars()); }

    private:
        std::vector<double> invoke(const std::vector<IntegerGene>& x) const override
        {
            assert(x.size() == num_vars());

            double fitness = 0.0;
            for (size_t i = 0; i < x.size(); i++)
            {
                fitness += double(char(x[i] + 32) == target_[i]);
            }

            return { fitness };
        }

        std::string target_;
    };

} // namespace genetic_algorithm::benchmark

#endif // !GA_BENCHMARK_INTEGER_HPP