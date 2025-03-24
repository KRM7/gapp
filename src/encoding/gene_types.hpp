/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_ENCODING_GENE_TYPES_HPP
#define GAPP_ENCODING_GENE_TYPES_HPP

#include "../utility/type_traits.hpp"
#include "../utility/type_list.hpp"
#include <cstdint>
#include <cstddef>

namespace gapp
{
    /** The gene type used in the binary-encoded genetic algorithm. @see BinaryGA */
    using BinaryGene = std::uint8_t;

    /** The gene type used in the real-encoded genetic algorithm. @see RCGA */
    using RealGene = double;

    /** The gene type used in the permutation-encoded genetic algorithm. @see PermutationGA */
    using PermutationGene = std::size_t;

    /** The gene type used in the integer-encoded genetic algorithm. @see IntegerGA */
    using IntegerGene = std::int64_t;

    /** The mixed gene type used for the mixed encodings. Made up of a number of unique component genes. @see MixedGA */
    template<typename... Ts>
    struct MixedGene
    {
        static_assert(sizeof...(Ts) >= 2, "A mixed gene must have at least 2 component gene types.");
        static_assert(detail::unique_types_v<Ts...>, "The component gene types of the mixed gene must be unique.");
        static_assert((!detail::is_specialization_of_v<Ts, MixedGene> && ...), "The component gene types can't also be mixed genes.");
    };


    /**
    * Type trait to check if a particular gene type is bounded or not.
    * This struct should be specialized for custom gene types if they are bounded,
    * otherwise they will be treated as if they were not bounded.
    *
    * Example:
    * ```
    * template<> struct is_bounded_gene<SomeBoundedGeneType> : std::true_type {};
    * ```
    * 
    * @see is_bounded_gene_v
    */
    template<typename GeneType>
    struct is_bounded_gene : std::false_type {};

    template<> struct is_bounded_gene<RealGene> : std::true_type {};
    template<> struct is_bounded_gene<IntegerGene> : std::true_type {};

    /**
    * Type trait to check if a particular gene type is bounded or not.
    * 
    * When defining new bounded gene types, the is_bounded_gene struct should be
    * specialized instead of specializing this variable template directly.
    * 
    * @see is_bounded_gene
    */
    template<typename GeneType>
    inline constexpr bool is_bounded_gene_v = is_bounded_gene<GeneType>::value;


    /**
    * Type trait to check if the gene type has at least 1 bounded component gene.
    * For simple gene types, this is equivalent to is_bounded_gene, for mixed gene
    * types, it is true if at least 1 of the component genes is bounded.
    * 
    * @see is_partially_bounded_gene_v, is_bounded_gene
    */
    template<typename GeneType>
    struct is_partially_bounded_gene : is_bounded_gene<GeneType> {};

    template<template<typename...> class MixedGene, typename... Ts>
    struct is_partially_bounded_gene<MixedGene<Ts...>> : std::disjunction<is_bounded_gene<Ts>...> {};

    /**
    * Type trait to check if the gene type has at least 1 bounded component gene.
    * For simple gene types, this is equivalent to is_bounded_gene_v, for mixed gene
    * types, it is true if at least 1 of the component genes is bounded.
    *
    * @see is_partially_bounded_gene, is_bounded_gene_v
    */
    template<typename GeneType>
    inline constexpr bool is_partially_bounded_gene_v = is_partially_bounded_gene<GeneType>::value;


    /** Type trait to check if a particular gene type is mixed or not. @see is_mixed_gene_v */
    template<typename GeneType>
    struct is_mixed_gene : detail::is_specialization_of<GeneType, MixedGene> {};

    /** Type trait to check if a particular gene type is mixed or not. @see is_mixed_gene */
    template<typename GeneType>
    inline constexpr bool is_mixed_gene_v = is_mixed_gene<GeneType>::value;


    /** A type list containing the component genes of GeneType. */
    template<typename GeneType>
    using component_genes_t = std::conditional_t<is_mixed_gene_v<GeneType>, detail::args_to_list_t<GeneType>, detail::type_list<GeneType>>;

    /** A type list containing the bounded component genes of GeneType. */
    template<typename GeneType>
    using bounded_component_genes_t = typename component_genes_t<GeneType>::template filter_types_t<is_bounded_gene>;

} // namespace gapp

#endif // !GAPP_ENCODING_GENE_TYPES_HPP
