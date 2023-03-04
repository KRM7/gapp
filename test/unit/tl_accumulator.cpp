/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/tl_accumulator.hpp"
#include "utility/iterators.hpp"
#include "utility/functional.hpp"
#include <vector>
#include <algorithm>
#include <execution>

using namespace genetic_algorithm::detail;

TEST_CASE("tl_vector_accumulator", "[tl_accumulator]")
{
    constexpr size_t nrows = 100'000;
    constexpr size_t ncols = 1'000;

    const std::vector mat(nrows, std::vector(ncols, 1));

    tl_vector_accumulator<int>::reset(ncols);

    std::for_each(std::execution::par_unseq, mat.begin(), mat.end(), [](const std::vector<int>& row)
    {
        for (size_t col = 0; col < row.size(); col++)
        {
            tl_vector_accumulator<int>::at(col) += row[col];
        }
    });

    const auto colwise_sums = tl_vector_accumulator<int>::collect();

    REQUIRE(
        std::all_of(colwise_sums.begin(), colwise_sums.end(), equal_to(int(nrows)))
    );
}