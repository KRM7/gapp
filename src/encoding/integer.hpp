/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_INTEGER_HPP
#define GA_ENCODING_INTEGER_HPP

#include "gene_types.hpp"
#include "../core/ga_base.hpp"
#include "../population/candidate.hpp"
#include "../crossover/integer.hpp"
#include "../mutation/integer.hpp"
#include <cstddef>

namespace genetic_algorithm
{
    template<>
    struct GaTraits<IntegerGene>
    {
        using DefaultCrossover = crossover::integer::TwoPoint;
        using DefaultMutation = mutation::integer::Uniform;

        static constexpr Probability defaultMutationRate(size_t chrom_len) noexcept { return 1.0 / chrom_len; }
    };

    /**
    * Integer encoded genetic algorithm. \n
    * Similar to the @ref BinaryGA, but the values of the genes can be any integer in a closed interval, not just 0 or 1.
    * The interval is specified by the gene bounds. \n
    * The algorithm also uses a modified mutation function with swaps and inversions.
    */
    class IntegerGA final : public GA<IntegerGene>
    {
    public:
        using GA::GA;
    private:
        Candidate<GeneType> generateCandidate() const override;
    };

} // namespace genetic_algorithm

#endif // !GA_ENCODING_INTEGER_HPP