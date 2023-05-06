/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_TRAITS_HPP
#define GA_CORE_GA_TRAITS_HPP

#include "../utility/type_traits.hpp"

namespace genetic_algorithm
{
    /**
    * Traits class describing some attributes of each GA type. \n
    * When defining a new class that inherits from the GA<T> class, this class must be specialized for the gene type T used
    * before the derived class is declared, and must have the following members:
    *   - DefaultCrossover type
    *   - DefaultMutation type
    *   - defaultMutationRate static member function
    * 
    * The following gene types are reserved for the GAs already implemented in the library
    * and can't be used as the gene type of new encodings:
    *   int, std::int8_t, std::size_t, double
    */
    template<typename GeneType>
    struct GaTraits
    {
        static_assert(detail::always_false<GeneType>, "The GaTraits class must be specialized for each new gene type.");

        // using DefaultCrossover = ...
        // using DefaultMutation = ...

        // static Probability defaultMutationRate(size_t chromosome_size) { ... }
    };

} // namespace genetic_algorithm

#endif // !GA_CORE_GA_TRAITS_HPP