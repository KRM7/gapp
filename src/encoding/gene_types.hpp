/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_GENE_TYPES_HPP
#define GA_ENCODING_GENE_TYPES_HPP

#include <cstdint>
#include <cstddef>

namespace gapp
{
    /** The gene type used in the binary-encoded genetic algorithm. @see BinaryGA */
    using BinaryGene = std::int8_t;

    /** The gene type used in the real-encoded genetic algorithm. @see RCGA */
    using RealGene = double;

    /** The gene type used in the permutation-encoded genetic algorithm. @see PermutationGA */
    using PermutationGene = std::size_t;

    /** The gene type used in the integer-encoded genetic algorithm. @see IntegerGA */
    using IntegerGene = std::int64_t;

    /**
    * Type trait to check if a particular gene type is bounded or not.
    * This variable should be specialized for custom gene types if they are bounded,
    * otherwise they will be treated as if they were not bounded.
    * 
    * Example:
    * ```
    * template<> inline constexpr bool is_bounded<MyGene> = true;
    * ```
    */
    template<typename GeneType>
    inline constexpr bool is_bounded = false;

    template<> inline constexpr bool is_bounded<RealGene>    = true;
    template<> inline constexpr bool is_bounded<IntegerGene> = true;

} // namespace gapp

#endif // !GA_ENCODING_GENE_TYPES_HPP