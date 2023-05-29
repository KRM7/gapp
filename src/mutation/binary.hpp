/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BINARY_HPP
#define GA_MUTATION_BINARY_HPP

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
    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;
    };

} // namespace gapp::mutation::binary

#endif // !GA_MUTATION_BINARY_HPP