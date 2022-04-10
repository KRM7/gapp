/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_INTEGER_HPP
#define GA_MUTATION_INTEGER_HPP

#include "mutation_base.hpp"

#include <cstddef>

/** Predefined mutation operators for the integer encoded genetic algorithm. */
namespace genetic_algorithm::mutation::integer
{
    /**
    * Uniform mutation operator for the integer encoded genetic algorithm.
    * Each gene of the chromosome is changed with the specified mutation probability
    * to another value selected from a uniform distribution over all other values.
    */
    class Uniform : public Mutation<size_t>
    {
    public:
        using Mutation::Mutation;
    private:
        void mutate(const GaInfo& ga, Candidate<size_t>& candidate) const override;
    };

} // namespace genetic_algorithm::mutation::integer

#endif // !GA_MUTATION_INTEGER_HPP