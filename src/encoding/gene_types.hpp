/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_ENCODING_GENE_TYPES_HPP
#define GA_ENCODING_GENE_TYPES_HPP

#include <cstdint>
#include <cstddef>

namespace genetic_algorithm
{
    using BinaryGene      = std::int8_t;
    using RealGene        = double;
    using PermutationGene = std::size_t;
    using IntegerGene     = std::int64_t;

} // namespace genetic_algorithm

#endif // !GA_ENCODING_GENE_TYPES_HPP