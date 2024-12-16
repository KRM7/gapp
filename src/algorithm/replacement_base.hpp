/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ALGORITHM_SOGA_REPLACEMENT_BASE_HPP
#define GAPP_ALGORITHM_SOGA_REPLACEMENT_BASE_HPP

#include "../core/population.hpp"
#include "../utility/small_vector.hpp"
#include <cstddef>

namespace gapp
{
    class GaInfo;

} // namespace gapp

/** Population replacement policies for the single-objective algorithms. */
namespace gapp::replacement
{
    /**
    * This is the base class used for all of the single-objective population replacement policies.
    * The replacement operator is used to select the candidates of the next population from the
    * combined parent and child populations.
    * 
    * New replacement policies for the single-objective algorithm should be derived from this
    * class, and there is a single method that should be implemented in the derived classes:
    * 
    *   - nextPopulationImpl :  Selects the candidates for the next population.
    */
    class Replacement
    {
    public:
        /**
        * Select the candidates of the next generation from the candidates of the
        * combined current and child populations.
        * 
        * The @pop parameter describes the combined parent and child populations.
        * The top half (first population_size elements) corresponds to the parent
        * population, while the rest (another population_size elements) corresponds
        * to the child population.
        * 
        * The method should return population_size number of unique indices into
        * the population @p pop. The candidates at these indices will comprise the
        * next generation's population.
        * 
        * @param ga The %GA that uses the replacement method.
        * @param pop A view of the combined parent and child population.
        * @returns The indices of the candidates selected from the population @p pop.
        */
        virtual small_vector<size_t> nextPopulationImpl(const GaInfo& ga, const PopulationView& pop) = 0;
        

        /** Destructor. */
        virtual ~Replacement()                      = default;

    protected:

        Replacement()                               = default;
        Replacement(const Replacement&)             = default;
        Replacement(Replacement&&)                  = default;
        Replacement& operator=(const Replacement&)  = default;
        Replacement& operator=(Replacement&&)       = default;
    };

} // namespace gapp::replacement

#endif // !GAPP_ALGORITHM_SOGA_REPLACEMENT_BASE_HPP
