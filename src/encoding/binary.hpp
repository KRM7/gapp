﻿/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_BINARY_HPP
#define GA_ENCODING_BINARY_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../core/ga_traits.hpp"
#include "../population/candidate.hpp"
#include "../crossover/binary.hpp"
#include "../mutation/binary.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    template<>
    struct GaTraits<BinaryGene>
    {
        using DefaultCrossover = crossover::binary::TwoPoint;
        using DefaultMutation = mutation::binary::Flip;

        static constexpr Probability defaultMutationRate(size_t chrom_len) noexcept { return 1.0 / chrom_len; }
    };

    /** Standard binary encoded genetic algorithm. */
    class BinaryGA final : public GA<BinaryGene>
    {
    public:
        using GA::GA;
    private:
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_BINARY_HPP