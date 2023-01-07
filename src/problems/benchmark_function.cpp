/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "benchmark_function.hpp"
#include "../core/ga_base.decl.hpp"
#include "../encoding/gene_types.hpp"
#include <vector>
#include <numeric>
#include <cmath>
#include <cassert>

namespace genetic_algorithm::problems
{
    using RealBoundsVec = std::vector<typename GA<RealGene>::GeneBounds>;

    static std::vector<RealGene> convert(const std::vector<BinaryGene>& binary_chrom, const RealBoundsVec& bounds, size_t var_bits)
    {
        std::vector<RealGene> vars(bounds.size());

        for (size_t i = 0; i < vars.size(); i++)
        {
            const auto first = binary_chrom.begin() + i * var_bits;
            const auto last = binary_chrom.begin() + (i + 1) * var_bits;

            const RealGene val = std::accumulate(first, last, 0.0, [](RealGene acc, BinaryGene bit) noexcept
            {
                return (acc * 2) + bit;
            });

            vars[i] = val / (std::pow(2.0, var_bits) - 1); // use double to avoid integer overflow
            vars[i] *= bounds[i].upper - bounds[i].lower;
            vars[i] += bounds[i].lower;
        }

        return vars;
    }

    std::vector<double> BenchmarkFunctionReal1::operator()(const std::vector<BinaryGene>& binary_chrom) const
    {
        assert((binary_chrom.size() / var_bits_) == num_vars());

        return convert(binary_chrom, bounds_, var_bits_);
    }

    std::vector<double> BenchmarkFunctionRealN::operator()(const std::vector<BinaryGene>& binary_chrom) const
    {
        assert((binary_chrom.size() / var_bits_) == num_vars());

        return convert(binary_chrom, bounds_, var_bits_);
    }

} // namespace genetic_algorithm::problems