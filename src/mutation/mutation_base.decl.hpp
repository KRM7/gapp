/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_MUTATION_BASE_DECL_HPP
#define GA_MUTATION_BASE_DECL_HPP

#include "../population/candidate.hpp"
#include "../utility/bounded_value.hpp"
#include <vector>
#include <utility>

namespace gapp
{
    template<typename T>
    class GA;

} // namespace gapp

namespace gapp::mutation
{
    /**
    * The base class used for the mutation operators of the GAs.
    * %Mutation operators take a candidate solution, and modify it in some way with a
    * given probability. This probability can be interpreted either per-candidate or
    * per-gene depending on how the particular operator is defined.
    * 
    * New mutation operators should be derived from this class, and they should
    * implement the following virtual methods:
    * 
    *   - mutate : Perform the mutation on a single candidate's chromosome.
    * 
    * @tparam T The gene type the mutation operator is defined for.
    */
    template<typename T>
    class Mutation
    {
    public:
        /** The gene type the mutation operator is defined for. */
        using GeneType = T;

        /**
        * Create a mutation operator.
        *
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Mutation(Probability pm) noexcept :
            pm_(pm)
        {}

        /**
        * Set the mutation rate used by the operator.
        * 
        * @param pm The mutation probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void mutation_rate(Probability pm) noexcept { pm_ = pm; };

        /** @returns The mutation rate set for the operator. */
        [[nodiscard]]
        constexpr Probability mutation_rate() const noexcept { return pm_; }

        /**
        * Perform mutation on a candidate using the set mutation probability.
        * Implemented by mutate().
        *
        * @param ga The genetic algorithm the mutation operator is being used in.
        * @param candidate The candidate to mutate.
        */
        void operator()(const GA<T>& ga, Candidate<T>& candidate) const;


        /** Destructor. */
        virtual ~Mutation()                     = default;

    protected:

        Mutation(const Mutation&)               = default;
        Mutation(Mutation&&)                    = default;
        Mutation& operator=(const Mutation&)    = default;
        Mutation& operator=(Mutation&&)         = default;

    private:

        /**
        * The implementation of the mutation operator. Performs the mutation
        * on the given chromosome in-place with the set probability.
        * This function must handle the mutation probability properly as part
        * of its implementation. The mutated chromosome should be valid candidate
        * solution for the given problem and %GA.
        * 
        * This method will be called exactly once for each child solution
        * in every population.
        *
        * The function must be thread-safe if parallel execution is enabled for the
        * GAs (which is true by default).
        * 
        * @param ga The genetic algorithm the mutation operator is being used in.
        * @param candidate The candidate solution that will be mutated.
        * @param chromosome The chromosome to mutate. This is the chromosome of @p candidate.
        */
        virtual void mutate(const GA<T>& ga, const Candidate<T>& candidate, Chromosome<T>& chromosome) const = 0;

        Probability pm_;
    };

} // namespace gapp::mutation

#endif // !GA_MUTATION_BASE_DECL_HPP