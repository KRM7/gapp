/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_CORE_GA_TRAITS_HPP
#define GA_CORE_GA_TRAITS_HPP

#include "../utility/type_traits.hpp"

namespace gapp
{
    /**
    * Traits class describing some attributes of each %GA class. \n
    * When defining a new genetic algorithm class that inherits from the %GA class,
    * this class must be specialized for the gene type T used before the derived class is
    * declared, and must have the following members:
    * 
    *   - DefaultCrossover type (must be default constructible)
    *   - DefaultMutation type  (must be constructible using the return value of defaultMutationRate)
    *   - defaultMutationRate(size_t chrom_len) -> Probability (static member function)
    * 
    * The following gene types are reserved for the GAs already implemented in the library
    * and can't be used as the gene type of new encodings:
    *   - std::int8_t
    *   - std::size_t
    *   - int
    *   - double
    * 
    * Example:
    * ```
    * template<>
    * struct GaTraits<MyGeneType>
    * {
    *   using DefaultCrossover = MyCrossoverType;
    *   using DefaultMutation  = MyMutationType;
    *   static Probability defaultMutationRate(size_t chromosome_size) { return 0.01; }
    * };
    * ```
    * 
    * @see GA
    */
    template<typename GeneType>
    struct GaTraits
    {
        static_assert(detail::always_false<GeneType>, "The GaTraits class must be specialized for each new gene type.");
    };

} // namespace gapp

#endif // !GA_CORE_GA_TRAITS_HPP