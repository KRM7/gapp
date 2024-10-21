/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include "benchmark_function.hpp"
#include "../encoding/gene_types.hpp"
#include "../utility/utility.hpp"
#include <numeric>
#include <cstddef>

namespace gapp::problems
{
    Candidate<RealGene> BenchmarkFunction<RealGene>::convert(const Candidate<BinaryGene>& sol, const BoundsVector<RealGene>& bounds, size_t var_bits) const
    {
        GAPP_ASSERT(sol.chromosome.size() == var_bits * bounds.size());

        Candidate<RealGene> new_sol(bounds.size());

        for (size_t i = 0; i < new_sol.size(); i++)
        {
            const auto first = sol.begin() + i * var_bits;
            const auto last = sol.begin() + (i + 1) * var_bits;

            const RealGene val = std::accumulate(first, last, 0.0, [&](RealGene acc, BinaryGene bit) noexcept
            {
                return (acc * 2) + bit * lsb_;
            });

            new_sol[i] = val * (bounds[i].upper() - bounds[i].lower()) + bounds[i].lower();
        }

        return new_sol;
    }

} // namespace gapp::problems
