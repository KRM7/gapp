/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_PERMUTATION_HPP
#define GA_ENCODING_PERMUTATION_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/permutation.hpp"
#include "../mutation/permutation.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    template<>
    struct GaTraits<PermutationGene>
    {
        using DefaultCrossover = crossover::perm::Order2;
        using DefaultMutation = mutation::perm::Inversion;

        static constexpr Probability defaultMutationRate(size_t) noexcept { return 0.6; }
    };

    /**
    * Genetic algorithm in which the chromosomes encode permutations. \n
    * The genes of the chromosomes are unique unsigned integers on the closed interval [0, chrom_len - 1].
    * 
    * The first and last elements of the permutations are assumed to be unrelated, eg. the permutation
    * A-B-C-D will not be considered equal to the permutation B-C-D-A by the algorithm.
    */
    class PermutationGA final : public GA<PermutationGene>
    {
    public:
        using GA::GA;
    private:
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_PERMUTATION_HPP