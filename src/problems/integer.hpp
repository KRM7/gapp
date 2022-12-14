/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/** Test problems for the integer-encoded algorithms. */

#ifndef GA_PROBLEMS_INTEGER_HPP
#define GA_PROBLEMS_INTEGER_HPP

#include "benchmark_function.hpp"
#include <vector>
#include <string>
#include <utility>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
{
    /**
    * Implementation of a simple test problem for the integer-encoded algorithm. \n
    * The goal is for the algorithm to match a target string.
    * 
    * The problem is implemented for maximization, and only usable with the single-objective,
    * integer-encoded algorithm. \n
    * The number of variables will be equal to the length of the target string.
    */
    class StringFinder : public BenchmarkFunction<IntegerGene>
    {
    public:
        /**
        * Construct an instance of the string matching problem.
        * 
        * @param target The string to look for.
        */
        explicit StringFinder(std::string target) :
            BenchmarkFunction<IntegerGene>("StringFinder", 1, target.size(), Bounds{ 0, 95 }), target_(std::move(target))
        {}

        /** @returns The fitness value of the optimal solution. */
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

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_INTEGER_HPP