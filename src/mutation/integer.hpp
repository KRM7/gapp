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
        /**
        * Create a uniform crossover operator with the specified parameters.
        * 
        * @param base The number of values a gene can take. Genes can be integers in the range [0, base-1].
        * @param The mutation probability used.
        */
        explicit Uniform(size_t base, double pm = 0.01);

        /**
        * Sets the base for this crossover operator.
        * 
        * @param base The number of values a gene can take. Genes can be integers in the range [0, base-1].
        */
        void base(size_t base);

        /** @returns The base parameter currently set for this operator. */
        [[nodiscard]]
        size_t base() const noexcept { return base_; }

    private:
        void mutate(const GaInfo& ga, Candidate<size_t>& candidate) const override;

        size_t base_;
    };

} // namespace genetic_algorithm::mutation::integer

#endif // !GA_MUTATION_INTEGER_HPP