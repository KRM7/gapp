/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/utility.hpp"
#include <numeric>
#include <cmath>
#include <cstddef>

namespace gapp::problems
{
    Chromosome<RealGene> BenchmarkFunction<RealGene>::convert(const Chromosome<BinaryGene>& bchrom, const BoundsVector<RealGene>& bounds, size_t var_bits)
    {
        GAPP_ASSERT(bchrom.size() == var_bits * bounds.size());

        Chromosome<RealGene> vars(bounds.size());

        for (size_t i = 0; i < vars.size(); i++)
        {
            const auto first = bchrom.begin() + i * var_bits;
            const auto last = bchrom.begin() + (i + 1) * var_bits;

            const RealGene val = std::accumulate(first, last, 0.0, [](RealGene acc, BinaryGene bit) noexcept
            {
                return (acc * 2) + bit;
            });

            vars[i] = val / (std::pow(2.0, var_bits) - 1); // use double to avoid integer overflow
            vars[i] *= bounds[i].upper() - bounds[i].lower();
            vars[i] += bounds[i].lower();
        }

        return vars;
    }

} // namespace gapp::problems