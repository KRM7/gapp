/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_DECL_HPP
#define GA_MUTATION_BASE_DECL_HPP

#include "../population/candidate.hpp"
#include <vector>
#include <utility>

namespace genetic_algorithm
{
    class GaInfo;

} // namespace genetic_algorithm

/** Mutation operators used in the algorithms. */
namespace genetic_algorithm::mutation
{
    /**
    * Base class used for all of the mutation operators.
    * Every implemented mutation operator takes a Candidate, and mutates that candidate using the set probability.
    * (This probability can either be per-candidate or per-gene depending on the operator.)
    */
    template<Gene T>
    class Mutation
    {
    public:
        using GeneType = T;

        /**
        * Create a mutation operator.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        explicit Mutation(double pm);

        Mutation(const Mutation&)               = default;
        Mutation(Mutation&&)                    = default;
        Mutation& operator=(const Mutation&)    = default;
        Mutation& operator=(Mutation&&)         = default;
        virtual ~Mutation()                     = default;

        /**
        * Set the mutation rate used for the operator to @p pm.
        * 
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        void mutation_rate(double pm);

        /** @returns The mutation rate currently set for this mutation operator. */
        [[nodiscard]]
        double mutation_rate() const noexcept { return pm_; }

        /**
        * Perform the mutation operation on a candidate using the set mutation rate.
        *
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param candidate The candidate to mutate.
        */
        void operator()(const GaInfo& ga, Candidate<T>& candidate) const;

    protected:

        /* The actual mutation function. Performs the mutation using the set probability. */
        virtual void mutate(const GaInfo& ga, Candidate<T>& candidate) const = 0;

        double pm_;
    };

} // namespace genetic_algorithm::mutation

#endif // !GA_MUTATION_BASE_DECL_HPP