/* Copyright (c) 2025 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "encoding/gene_types.hpp"

using namespace gapp;

struct CustomGeneType {};
struct CustomBoundedGeneType {};

namespace gapp
{
    template<> struct is_bounded_gene<CustomBoundedGeneType> : std::true_type {};
}

TEST_CASE("is_bounded_gene", "[gene_types]")
{
    STATIC_REQUIRE(is_bounded_gene_v<RealGene>);
    STATIC_REQUIRE(!is_bounded_gene_v<BinaryGene>);
    STATIC_REQUIRE(!is_bounded_gene_v<PermutationGene>);

    STATIC_REQUIRE(is_bounded_gene_v<CustomBoundedGeneType>);
    STATIC_REQUIRE(!is_bounded_gene_v<CustomGeneType>);
}

TEST_CASE("is_partially_bounded_gene", "[gene_types]")
{
    STATIC_REQUIRE(is_partially_bounded_gene_v<RealGene>);
    STATIC_REQUIRE(!is_partially_bounded_gene_v<BinaryGene>);

    STATIC_REQUIRE(is_partially_bounded_gene_v<CustomBoundedGeneType>);
    STATIC_REQUIRE(!is_partially_bounded_gene_v<CustomGeneType>);

    STATIC_REQUIRE(is_partially_bounded_gene_v<MixedGene<RealGene, BinaryGene>>);
    STATIC_REQUIRE(!is_partially_bounded_gene_v<MixedGene<PermutationGene, BinaryGene>>);

    STATIC_REQUIRE(is_partially_bounded_gene_v<MixedGene<BinaryGene, PermutationGene, CustomBoundedGeneType>>);
    STATIC_REQUIRE(!is_partially_bounded_gene_v<MixedGene<BinaryGene, PermutationGene, CustomGeneType>>);
}

TEST_CASE("is_mixed_gene", "[gene_types]")
{
    STATIC_REQUIRE(!is_mixed_gene_v<RealGene>);
    STATIC_REQUIRE(!is_mixed_gene_v<BinaryGene>);
    STATIC_REQUIRE(!is_mixed_gene_v<CustomGeneType>);

    STATIC_REQUIRE(is_mixed_gene_v<MixedGene<RealGene, BinaryGene>>);
    STATIC_REQUIRE(is_mixed_gene_v<MixedGene<CustomGeneType, BinaryGene, CustomBoundedGeneType>>);
}

TEST_CASE("component_genes_t", "[gene_types]")
{
    STATIC_REQUIRE(std::is_same_v<component_genes_t<RealGene>, detail::type_list<RealGene>>);
    STATIC_REQUIRE(std::is_same_v<component_genes_t<BinaryGene>, detail::type_list<BinaryGene>>);

    STATIC_REQUIRE(std::is_same_v<component_genes_t<MixedGene<RealGene, BinaryGene>>, detail::type_list<RealGene, BinaryGene>>);
}

TEST_CASE("bounded_component_genes_t", "[gene_types]")
{
    STATIC_REQUIRE(std::is_same_v<bounded_component_genes_t<RealGene>, detail::type_list<RealGene>>);
    STATIC_REQUIRE(std::is_same_v<bounded_component_genes_t<BinaryGene>, detail::type_list<>>);

    STATIC_REQUIRE(std::is_same_v<bounded_component_genes_t<MixedGene<RealGene, BinaryGene>>, detail::type_list<RealGene>>);
}
