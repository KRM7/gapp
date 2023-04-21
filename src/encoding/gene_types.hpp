/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_GENE_TYPES_HPP
#define GA_ENCODING_GENE_TYPES_HPP

#include <cstdint>
#include <cstddef>

namespace genetic_algorithm
{
    using BinaryGene      = std::int8_t;    /**< The gene type used in the binary encoded algorithms (BinaryGA). */
    using RealGene        = double;         /**< The gene type used in the real encoded algorithms (RCGA). */
    using PermutationGene = std::size_t;    /**< The gene type used in the permutation encoded algorithms (PermutationGA). */
    using IntegerGene     = std::int64_t;   /**< The gene type used in the integer encoded algorithms (IntegerGA). */

    /**
    * Type trait to check if a particular gene type is bounded or not. \n
    * Must be specialized for custom gene types if they are bounded.
    */
    template<typename GeneType>
    inline constexpr bool is_bounded = false;

    template<> inline constexpr bool is_bounded<RealGene>    = true;
    template<> inline constexpr bool is_bounded<IntegerGene> = true;

} // namespace genetic_algorithm

#endif // !GA_ENCODING_GENE_TYPES_HPP