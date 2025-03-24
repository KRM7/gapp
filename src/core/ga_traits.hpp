/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_CORE_GA_TRAITS_HPP
#define GAPP_CORE_GA_TRAITS_HPP

#include "../utility/type_traits.hpp"

namespace gapp
{
    /**
    * Traits class describing some attributes of each %GA class. \n
    * When defining a new genetic algorithm class that inherits from the %GA class,
    * this class must be specialized for the gene type T used, before the derived class
    * is declared, and must have the following members:
    * 
    *   - DefaultCrossover type (must be default constructible)
    *   - DefaultMutation type  (must be constructible using the return value of defaultMutationRate)
    *   - defaultMutationRate(size_t chrom_len) -> Probability (static member function)
    * 
    * If T is not a bounded gene type:
    *   - randomChromosome(size_t chrom_len) -> Chromosome<T> (static member function)
    * 
    * or, if T is a bounded gene type:
    *   - randomChromosome(size_t chrom_len, BoundsView<T>) -> Chromosome<T> (static member function)
    * 
    * The following gene types are reserved for the GAs already implemented in the library
    * and can't be used as the gene type of new encodings:
    *   - std::uint8_t
    *   - std::size_t
    *   - std::int64_t
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
    * @note This class shouldn't be specialized for individual instances of the MixedGene
    *   gene type, but a specialization must exist for each component gene of the mixed gene
    *   (i.e. for each type parameter of the MixedGene instance). This means that if a mixed
    *   gene uses a custom gene type as one of it's component genes, the class should only
    *   be specialized for that custom gene type.
    * 
    * @see GA
    */
    template<typename GeneType>
    struct GaTraits
    {
        static_assert(detail::always_false<GeneType>, "The GaTraits class must be specialized for every gene type.");
    };

} // namespace gapp

#endif // !GAPP_CORE_GA_TRAITS_HPP
