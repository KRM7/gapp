/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BINARY_HPP
#define GA_MUTATION_BINARY_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"

/** Predefined mutation operators for the binary encoded genetic algorithm. */
namespace genetic_algorithm::mutation::binary
{
    /**
    * Standard flip mutation for the binary encoded genetic algorithm.
    * Each gene of the chromosome is flipped (0 -> 1, or 1 -> 0) with the specified mutation probability.
    */
    class Flip final : public Mutation<BinaryGene>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GaInfo& ga, Candidate<GeneType>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::binary

#endif // !GA_MUTATION_BINARY_HPP