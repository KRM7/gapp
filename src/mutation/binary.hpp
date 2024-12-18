/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_MUTATION_BINARY_HPP
#define GAPP_MUTATION_BINARY_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"

/** Predefined mutation operators for the binary encoded genetic algorithm. */
namespace gapp::mutation::binary
{
    /**
    * Standard flip mutation for the binary encoded genetic algorithm.
    * Each gene of the chromosome is flipped (either from 0 to 1, or from 1 to 0)
    * with the specified mutation probability.
    */
    class Flip final : public Mutation<BinaryGene>
    {
    public:
        using Mutation::Mutation;
        constexpr bool allow_variable_chrom_length() const noexcept override { return true; }
    private:
        void mutate(const GaInfo& ga, const Candidate<GeneType>& candidate, Chromosome<GeneType>& chromosome) const override;
    };

} // namespace gapp::mutation::binary

#endif // !GAPP_MUTATION_BINARY_HPP
