/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CROSSOVER_BASE_DECL_HPP
#define GA_CROSSOVER_BASE_DECL_HPP

#include "../core/candidate.hpp"
#include "../utility/bounded_value.hpp"
#include <vector>
#include <utility>

namespace gapp
{
    template<typename T>
    class GA;

} // namespace gapp

namespace gapp::crossover
{
    /**
    * The base class used for the crossover operators of the GAs.
    * 
    * %Crossover operators take 2 candidate solutions (the parents), and create 2 new
    * candidates (children) based on the the parent candidates.
    * The crossover operation is only performed on the 2 parents with a set probability only,
    * the rest of the time the returned children will be the same as the parents.
    * 
    * New crossover operators should be derived from this class, and they must implement the
    * following virtual methods:
    * 
    *   - crossover : Perform the crossover on 2 candidate solutions.
    * 
    * @tparam The gene type the crossover operator is defined for.
    */
    template<typename T>
    class Crossover
    {
    public:
        /** The gene type the crossover operator is defined for. */
        using GeneType = T;

        /** Create a crossover operator using the default crossover probability. */
        constexpr explicit Crossover() noexcept :
            pc_(0.8_p) {}

        /**
        * Create a crossover operator.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr explicit Crossover(Probability pc) noexcept :
            pc_(pc) {}

        /**
        * Set the crossover probability used for the crossovers.
        *
        * @param pc The crossover probability. Must be in the closed interval [0.0, 1.0].
        */
        constexpr void crossover_rate(Probability pc) noexcept { pc_ = pc; }

        /** @returns The crossover rate set for the operator. */
        [[nodiscard]]
        constexpr Probability crossover_rate() const noexcept { return pc_; }

        /**
        * Perform the crossover operation on 2 candidate solutions with the set probability.
        * This function is implemented by crossover().
        *
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param parent1 The first parent solution.
        * @param parent2 The second parent solution.
        * @returns The pair of children resulting from the crossover.
        */
        CandidatePair<T> operator()(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const;


        /** Destructor. */
        virtual ~Crossover()                    = default;

    protected:

        Crossover(const Crossover&)             = default;
        Crossover(Crossover&&)                  = default;
        Crossover& operator=(const Crossover&)  = default;
        Crossover& operator=(Crossover&&)       = default;

    private:

        /**
        * The implementation of the crossover operator. Performs the crossover operation
        * on 2 parent solutions to generate 2 child solutions from them.
        * The implementation of this function shouldn't handle the crossover probability,
        * but instead it should perform the crossover unconditionally.
        * The chromosomes of the returned children should be valid solutions for the given
        * problem and %GA, but the rest of their properties (eg. fitness) are irrelevant.
        * 
        * This method will be called once for every 2 children that need to be generated
        * (ie. population_size/2 number of times, rounded up if the population size is odd)
        * in every generation.
        * 
        * The function must be thread-safe if parallel execution is enabled for the
        * GAs (which is true by default).
        * 
        * @param ga The genetic algorithm the crossover operator is being used in.
        * @param parent1 The first parent solution.
        * @param parent2 The second parent solution.
        * @returns The pair of children resulting from the crossover.
        */
        virtual CandidatePair<T> crossover(const GA<T>& ga, const Candidate<T>& parent1, const Candidate<T>& parent2) const = 0;

        Probability pc_;
    };

} // namespace gapp::crossover

#endif // !GA_CROSSOVER_BASE_DECL_HPP