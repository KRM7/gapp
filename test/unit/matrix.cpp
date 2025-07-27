/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#include <catch2/catch_test_macros.hpp>
#include "utility/matrix.hpp"
#include <type_traits>
#include <algorithm>

using namespace gapp::detail;

TEST_CASE("matrix", "[matrix]")
{
    SECTION("row types")
    {
        using Row = Matrix<double>::RowRef;
        using ConstRow = Matrix<double>::ConstRowRef;

        STATIC_REQUIRE(std::is_same_v<Row::value_type, double>);
        STATIC_REQUIRE(std::is_same_v<Row::reference, double&>);
        STATIC_REQUIRE(std::is_same_v<Row::pointer, double*>);

        STATIC_REQUIRE(std::is_same_v<ConstRow::value_type, double>);
        STATIC_REQUIRE(std::is_same_v<ConstRow::reference, const double&>);
        STATIC_REQUIRE(std::is_same_v<ConstRow::pointer, const double*>);
    }

    SECTION("iterator types")
    {
        using Iterator = Matrix<double>::iterator;
        using ConstIterator = Matrix<double>::const_iterator;

        STATIC_REQUIRE(std::is_same_v<Iterator::value_type, Matrix<double>::value_type>);
        STATIC_REQUIRE(std::is_same_v<Iterator::reference, Matrix<double>::reference>);

        STATIC_REQUIRE(std::is_same_v<ConstIterator::value_type, Matrix<double>::value_type>);
        STATIC_REQUIRE(std::is_same_v<ConstIterator::reference, Matrix<double>::const_reference>);
    }

    Matrix<int> mat1;
    Matrix<int> mat2 = { { 1, 2, 3 }, { 4, 5, 6 } };
    Matrix<int> mat3(4, 3, 1); 

    const Matrix<int> cmat = mat2;

    SECTION("member access operators")
    {
        REQUIRE(mat2(0, 0) == 1);
        REQUIRE(mat2(0, 2) == 3);
        REQUIRE(mat2(1, 1) == 5);

        REQUIRE(mat2[0][0] == 1);
        REQUIRE(mat2[0][2] == 3);
        REQUIRE(mat2[1][1] == 5);

        REQUIRE(cmat(0, 1) == 2);
        REQUIRE(cmat[0][1] == 2);
    }

    SECTION("sizes")
    {
        REQUIRE(mat1.empty());
        REQUIRE(mat1.size() == 0);
        REQUIRE(mat1.nrows() == 0);
        REQUIRE(mat1.ncols() == 0);

        REQUIRE(!mat2.empty());
        REQUIRE(mat2.size() == 2);
        REQUIRE(mat2.nrows() == 2);
        REQUIRE(mat2.ncols() == 3);

        REQUIRE(mat3.nrows() == 4);
        REQUIRE(mat3.ncols() == 3);
    }

    SECTION("resize")
    {
        mat2.resize(3, 2);

        REQUIRE(mat2.nrows() == 3);
        REQUIRE(mat2.ncols() == 2);

        REQUIRE(mat2[1][1] == 4);

        mat2.resize(1, 6);

        REQUIRE(mat2.nrows() == 1);
        REQUIRE(mat2.ncols() == 6);

        REQUIRE(mat2[0][5] == 6);
    }

    SECTION("matrix comparisons")
    {
        REQUIRE(mat2 != mat1);
        REQUIRE(mat2 != mat3);
        REQUIRE(mat2 == Matrix{ { 1, 2, 3 }, { 4, 5, 6 } });
        REQUIRE(mat2 == cmat);
    }

    SECTION("append rows")
    {
        const std::vector row = { 1, 1, 1 };

        // vectors
        mat2.append_row(row);

        REQUIRE(mat2.nrows() == 3);
        REQUIRE(mat2 == Matrix{ { 1, 2, 3 }, { 4, 5, 6 }, { 1, 1, 1 } });

        mat2.append_row(std::vector{ 2, 2, 2 });

        REQUIRE(mat2.nrows() == 4);
        REQUIRE(mat2 == Matrix{ { 1, 2, 3 }, { 4, 5, 6 }, { 1, 1, 1 }, { 2, 2, 2 } });

        // another row
        mat2.append_row(mat3[0]);
        
        REQUIRE(mat2.nrows() == 5);
        REQUIRE(mat2 == Matrix{ { 1, 2, 3 }, { 4, 5, 6 }, { 1, 1, 1 }, { 2, 2, 2 }, { 1, 1, 1 } });

        // append to empty matrix 
        mat1.append_row(std::vector{ 2, 3, 4, 1 });

        REQUIRE(mat1.ncols() == 4);
        REQUIRE(mat1.nrows() == 1);
        REQUIRE(mat1 == Matrix{ { 2, 3, 4, 1 } });
    }

    SECTION("swap")
    {
        using std::swap;
        swap(mat1, mat2);

        REQUIRE(mat2.empty());
        REQUIRE(mat1 == Matrix{ { 1, 2, 3 }, { 4, 5, 6 } });
    }
}

TEST_CASE("matrix_rows", "[matrix]")
{
    SECTION("row iterator types")
    {
        using Row = Matrix<double>::RowRef;
        using ConstRow = Matrix<double>::ConstRowRef;

        using Iterator = Row::iterator;
        using ConstIterator = ConstRow::iterator;

        STATIC_REQUIRE(std::is_same_v<Iterator::value_type, double>);
        STATIC_REQUIRE(std::is_same_v<Iterator::reference, double&>);
        STATIC_REQUIRE(std::is_same_v<Iterator::pointer, double*>);

        STATIC_REQUIRE(std::is_same_v<ConstIterator::value_type, double>);
        STATIC_REQUIRE(std::is_same_v<ConstIterator::reference, const double&>);
        STATIC_REQUIRE(std::is_same_v<ConstIterator::pointer, const double*>);
    }

    Matrix mat = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };

    const Matrix<int>& cmat = mat;

    const auto row1 = mat[0];
    const auto row2 = mat[1];
    const auto row3 = mat[2];
    const auto crow1 = cmat[0];
    const auto crow2 = cmat[1];
    const auto crow3 = cmat[2];

    SECTION("member access")
    {
        REQUIRE(row1[1] == 2);
        REQUIRE(row3[2] == 9);

        REQUIRE(crow2[0] == 4);
        REQUIRE(crow3[1] == 8);
    }

    SECTION("row comparisons")
    {
        REQUIRE(row1 == row1);
        REQUIRE(crow1 == row1);

        REQUIRE(row2 == std::vector{ 4, 5, 6 });
        REQUIRE(std::vector{ 7, 8, 9 } == crow3);
    }

    SECTION("sizes")
    {
        REQUIRE(row1.size() == 3);
        REQUIRE(crow2.size() == 3);

        REQUIRE(row3.ncols() == 3);
        REQUIRE(crow1.ncols() == 3);
    }

    SECTION("conversion")
    {
        Matrix<int>::ConstRowRef const_copy = row1;

        REQUIRE(const_copy == std::vector{ 1, 2, 3 });
    }

    SECTION("member assignment")
    {
        row1[1] = 13;
        REQUIRE(mat == Matrix{ { 1, 13, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });

        row2[0] = 31;
        REQUIRE(mat == Matrix{ { 1, 13, 3 }, { 31, 5, 6 }, { 7, 8, 9 } });
    }

    SECTION("row assignment")
    {
        row1 = crow2;
        REQUIRE(mat == Matrix{ { 4, 5, 6 }, { 4, 5, 6 }, { 7, 8, 9 } });

        row1 = row3;
        REQUIRE(cmat == Matrix{ { 7, 8, 9 }, { 4, 5, 6 }, { 7, 8, 9 } });

        row2 = std::vector{ 3, 2, 1 };
        REQUIRE(cmat == Matrix{ { 7, 8, 9 }, { 3, 2, 1 }, { 7, 8, 9 } });
    }

    SECTION("row swaps")
    {
        // swap rows
        using std::swap;
        swap(row1, row2);

        REQUIRE(mat == Matrix{ { 4, 5, 6 }, { 1, 2, 3 }, { 7, 8, 9 } });

        // swap row/vector
        std::vector vec = { 1, 2, 3 };
        swap(vec, row1);

        REQUIRE(mat == Matrix{ { 1, 2, 3 }, { 1, 2, 3 }, { 7, 8, 9 } });
        REQUIRE(vec == std::vector{ 4, 5, 6 });
    }

    SECTION("row swaps temp")
    {
        Matrix<int>::value_type temp = row1;
        row1 = row2;
        row2 = temp;

        REQUIRE(mat == Matrix{ { 4, 5, 6 }, { 1, 2, 3 }, { 7, 8, 9 } });
    }
}

TEST_CASE("matrix_algorithms", "[matrix]")
{
    Matrix mat1 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
    Matrix mat2 = { { 37, 40, 13 }, { 14, 4, 0 }, { 8, -1, 9 } };

    std::copy(mat2.begin(), mat2.end(), mat1.begin());
    REQUIRE(std::equal(mat1.begin(), mat1.end(), mat2.cbegin(), mat2.cend()));

    std::sort(mat1.begin(), mat1.end(), [](const auto& lhs, const auto& rhs) { return lhs[0] < rhs[0]; });
    REQUIRE(mat1 == Matrix{ { 8, -1, 9 }, { 14, 4, 0 }, { 37, 40, 13 } });

    std::reverse(mat1.begin(), mat1.end());
    REQUIRE(mat1 == Matrix{ { 37, 40, 13 }, { 14, 4, 0 }, { 8, -1, 9 } });
}
