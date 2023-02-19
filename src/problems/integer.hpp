/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

/** Test problems for the integer-encoded algorithms. */

#ifndef GA_PROBLEMS_INTEGER_HPP
#define GA_PROBLEMS_INTEGER_HPP

#include "benchmark_function.hpp"
#include <vector>
#include <string>
#include <cstddef>

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
    class StringFinder final : public BenchmarkFunction<IntegerGene>
    {
    public:
        /**
        * Create an instance of the string matching problem.
        * 
        * @param target The string to look for.
        */
        explicit StringFinder(std::string target);

    private:
        FitnessVector invoke(const std::vector<IntegerGene>& chrom) const override;

        std::string target_;
    };

} // namespace genetic_algorithm::problems

#endif // !GA_PROBLEMS_INTEGER_HPP