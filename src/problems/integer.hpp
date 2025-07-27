/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_PROBLEMS_INTEGER_HPP
#define GAPP_PROBLEMS_INTEGER_HPP

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include <string>

namespace gapp::problems
{
    /**
    * Implementation of a simple test problem for the integer-encoded %GA.
    * The goal is for the algorithm to match a target string.
    * 
    * The problem is implemented for maximization, and only usable with the single-objective,
    * integer-encoded %GA.
    * The number of variables will be equal to the length of the target string set in the
    * constructor.
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
        FitnessVector invoke(const Candidate<IntegerGene>& sol) const override;

        std::string target_;
    };

} // namespace gapp::problems

#endif // !GAPP_PROBLEMS_INTEGER_HPP
