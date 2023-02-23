/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#include "integer.hpp"
#include <vector>
#include <string>
#include <utility>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::problems
{
    StringFinder::StringFinder(std::string target) :
        BenchmarkFunction("StringFinder", target.size(), 1, Bounds{ 32, 32 + 95 }),
        target_(std::move(target))
    {
        optimal_value_ = { double(num_vars()) };
        ideal_point_ = optimal_value_;
        nadir_point_ = optimal_value_;
    }

    auto StringFinder::invoke(const std::vector<IntegerGene>& chrom) const -> FitnessVector
    {
        assert(chrom.size() == num_vars());

        double fitness = 0.0;
        for (size_t i = 0; i < chrom.size(); i++)
        {
            bool match = (static_cast<char>(chrom[i]) == target_[i]);
            fitness += double(match);
        }

        return { fitness };
    }

} // namespace genetic_algorithm::problems