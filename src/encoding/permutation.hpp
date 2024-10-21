/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ENCODING_PERMUTATION_HPP
#define GAPP_ENCODING_PERMUTATION_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../core/candidate.hpp"
#include "../crossover/permutation.hpp"
#include "../mutation/permutation.hpp"
#include <cstddef>

namespace gapp
{
    template<>
    struct GaTraits<PermutationGene>
    {
        using DefaultCrossover = crossover::perm::Order2;
        using DefaultMutation = mutation::perm::Inversion;

        static constexpr Probability defaultMutationRate(size_t) noexcept { return 0.6; }
    };

    /**
    * Permutation-encoded genetic algorithm class. This is the main solver
    * that should be used for combinatorial problems.
    * 
    * The chromosome of a candidate solution encodes a permutation. Every gene of a chromosome
    * is a unique unsigned integer in the closed interval of [0, chrom_len - 1].
    * 
    * Without any loss of generality, the first and last elements of the permutations are
    * assumed to be unrelated, e.g. the permutation A-B-C-D will not be considered equal
    * to the permutation B-C-D-A by the %GA. The fitness function should also be written
    * with this in mind.
    */
    class PermutationGA final : public GA<PermutationGene>
    {
    public:
        using GA::GA;
    private:
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace gapp

#endif // !GAPP_ENCODING_PERMUTATION_HPP
