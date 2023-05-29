/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_INTEGER_HPP
#define GA_MUTATION_INTEGER_HPP

#include "mutation_base.hpp"
#include "../encoding/gene_types.hpp"

/** Predefined mutation operators for the integer encoded genetic algorithm. */
namespace gapp::mutation::integer
{
    /**
    * %Uniform mutation operator for the integer encoded genetic algorithm.
    * Each gene of the chromosome is changed, with the specified mutation probability,
    * to another value selected from a uniform distribution over all other values.
    */
    class Uniform final : public Mutation<IntegerGene>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GA<GeneType>& ga, Candidate<GeneType>& candidate) const override;
    };

} // namespace gapp::mutation::integer

#endif // !GA_MUTATION_INTEGER_HPP