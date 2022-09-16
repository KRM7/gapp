/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATRIX_HPP
#define GA_UTILITY_MATRIX_HPP

#include "iterators.hpp"
#include <vector>
#include <memory>
#include <iterator>
#include <type_traits>
#include <utility>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<typename T, typename A = std::allocator<T>>
    class Matrix
    {
    public:
        using storage_type      = typename std::vector<T, A>;
        using value_type        = typename storage_type::value_type;
        using allocator_type    = typename storage_type::allocator_type;

        using size_type         = typename storage_type::size_type;
        using difference_type   = typename storage_type::difference_type;

        using reference         = typename storage_type::reference;
        using pointer           = typename storage_type::pointer;
        using const_reference   = typename storage_type::const_reference;
        using const_pointer     = typename storage_type::const_pointer;

        class Row;
        class RowIterator;
        class ConstRowIterator;

        using iterator               = RowIterator;
        using const_iterator         = ConstRowIterator;
        using reverse_iterator       = std::reverse_iterator<RowIterator>;
        using const_reverse_iterator = std::reverse_iterator<ConstRowIterator>;


        /* Special member functions */

        Matrix()                         = default;
        Matrix(const Matrix&)            = default;
        Matrix(Matrix&&)                 = default;
        Matrix& operator=(const Matrix&) = default;
        Matrix& operator=(Matrix&&)      = default;
        ~Matrix()                        = default;

        Matrix(const A& alloc) :
            nrows_(0), ncols_(0), data_(alloc)
        {}

        explicit Matrix(size_t square_size, const T& init = T(), const A& alloc = A()) :
            nrows_(square_size), ncols_(square_size), data_(square_size * square_size, init, alloc)
        {}

        Matrix(size_t nrows, size_t ncols, const T& init = T(), const A& alloc = A()) :
            nrows_(nrows), ncols_(ncols), data_(nrows * ncols, init, alloc)
        {}


        /* Operators */

        friend bool operator==(const Matrix&, const Matrix&) = default;

        Row operator[](size_t row) noexcept
        {
            return Row(*this, row);
        }

        const Row operator[](size_t row) const noexcept
        {
            return Row(*this, row);
        }


        /* Member access */

        reference operator()(size_type row, size_type col) noexcept
        {
            assert(row < nrows_ && col < ncols_);
            return data_[row * ncols_ + col];
        }

        const_reference operator()(size_type row, size_type col) const noexcept
        {
            assert(row < nrows_ && col < ncols_);
            return data_[row * ncols_ + col];
        }

        reference at(size_type row, size_type col)
        {
            return data_.at(row * ncols_ + col);
        }

        const_reference at(size_type row, size_type col) const
        {
            return data_.at(row * ncols_ + col);
        }

        pointer data() noexcept             { return data_.data(); }
        const_pointer data() const noexcept { return data_.data(); }


        /* Modifiers (for rows) */


        void push_back(const std::vector<T>& row)
        {
            assert(row.size() == ncols_);

            data_.reserve(data_.size() + row.size());
            data_.insert(data_.end(), row.begin(), row.end());
            nrows_++;
        }

        void push_back(std::vector<T>&& row)
        {
            assert(row.size() == ncols_);

            data_.reserve(data_.size() + row.size());
            data_.insert(data_.end(), std::make_move_iterator(row.begin()), std::make_move_iterator(row.end()));
            nrows_++;
        }

        void push_back(const Row& row)
        {
            assert(row.size() == ncols_);

            data_.reserve(data_.size() + row.size());
            data_.insert(data_.end(), row.begin(), row.end());
            nrows_++;
        }

        void pop_back()
        {
            assert(nrows_ != 0);

            data_.resize(--nrows_ * ncols_);
        }


        /* Size / capacity */

        size_type size() const noexcept { return nrows_; }
        size_type nrows() const noexcept { return nrows_; }
        size_type ncols() const noexcept { return ncols_; }
        size_type empty() const noexcept { return data_.empty(); }

        void reserve(size_type nrows, size_type ncols) { data_.reserve(nrows * ncols); }
        void resize(size_type nrows, size_type ncols) { data_.resize(nrows * ncols); nrows_ = nrows; ncols_ = ncols; }


        /* Other */

        void swap(Matrix& other) noexcept(std::is_nothrow_swappable_v<T>)
        {
            std::swap(data_, other.data_);
            std::swap(nrows_, other.nrows_);
            std::swap(ncols_, other.ncols_);
        }

    private:
        storage_type data_;
        size_type nrows_;
        size_type ncols_;
    };

    template<typename T, typename A>
    inline void swap(Matrix<T, A>& lhs, Matrix<T, A>& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }


    /* MATRIX ROW PROXY IMPLEMENTATION */

    template<typename T, typename A>
    class Matrix<T, A>::Row : public reverse_iterator_interface<Matrix<T, A>::Row>
    {
    private:
        Matrix& mat_;
        size_t row_;

    public:
        using value_type        = typename storage_type::value_type;
        using size_type         = typename storage_type::size_type;
        using difference_type   = typename storage_type::difference_type;
        using reference         = typename storage_type::reference;
        using const_reference   = typename storage_type::const_reference;


        Row(Matrix& mat, size_t row) noexcept :
            mat_(mat), row_(row) { assert(mat.nrows() > row); }

        Row(const Row&) = default;
        Row(Row&&)      = default;

        ~Row() = default;


        using iterator               = typename std::vector<T>::iterator;
        using const_iterator         = typename std::vector<T>::const_iterator;
        using reverse_iterator       = typename std::vector<T>::reverse_iterator;
        using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

        iterator begin() noexcept              { return mat_.begin() + row_ * mat_.ncols(); }
        const_iterator begin() const noexcept  { return mat_.begin() + row_ * mat_.ncols(); }
        iterator end() noexcept                { return mat_.begin() + (row_ + 1) * mat_.ncols(); }
        const_iterator end() const noexcept    { return mat_.begin() + (row_ + 1) * mat_.ncols(); }   


        Row& operator*() const noexcept { return *this; }

        reference operator[](size_t col) noexcept
        {
            assert(mat_.ncols() > col);
            return mat_(row_, col);
        }

        const_reference operator[](size_t col) const noexcept
        {
            assert(mat_.ncols() > col);
            return mat_(row_, col);
        }

        reference at(size_t col)
        {
            return mat_.at(row_, col);
        }

        const_reference at(size_t col) const
        {
            return mat_.at(row_, col);
        }

        size_type size() const noexcept
        {
            return mat_.ncols();
        }


        Row& operator=(const Row& rhs) noexcept(std::is_nothrow_assignable_v<T, T>)
        {
            assert(rhs.size() == this->size());

            for (size_t col = 0; col < rhs.size(); col++)
            {
                mat_(row_, col) = rhs[col];
            }

            return *this;
        }

        Row& operator=(const std::vector<T>& rhs) noexcept(std::is_nothrow_assignable_v<T, T>)
        {
            assert(rhs.size() == this->size());

            for (size_t col = 0; col < rhs.size(); col++)
            {
                mat_(row_, col) = rhs[col];
            }

            return *this;
        }

        /* We dont know if the underlying matrix can be moved from */
        Row& operator=(Row&& rhs) = delete;

        Row& operator=(std::vector<T>&& rhs)
        {
            assert(rhs.size() == this->size());

            for (size_t col = 0; col < rhs.size(); col++)
            {
                mat_(row_, col) = std::move_if_noexcept(rhs[col]);
            }

            return *this;
        }


        void swap(Row& other) noexcept(std::is_nothrow_swappable_v<T>)
        {
            assert(other.size() == this->size());

            for (size_t col = 0; col < other.size(); col++)
            {
                using std::swap;
                swap(mat_(row_, col), other[col]);
            }
        }

        void swap(std::vector<T>& other) noexcept(std::is_nothrow_swappable_v<T>)
        {
            assert(other.size() == this->size());

            for (size_t col = 0; col < other.size(); col++)
            {
                using std::swap;
                swap(mat_(row_, col), other[col]);
            }
        }

        explicit operator std::vector<T>() const
        {
            return std::vector<T>(this->begin(), this->end());
        }

        bool operator==(const Row& other) const noexcept(noexcept(std::declval<T>() != std::declval<T>()))
        {
            if ((&mat_ == &other.mat_) && (row_ == other.row_)) return true;

            if (this->size() != other.size()) return false;

            for (size_t col = 0; col < other.size(); col++)
            {
                if (mat_(row_, col) != other[col]) return false;
            }

            return true;
        }

        bool operator!=(const Row& other) const noexcept(noexcept(std::declval<T>() != std::declval<T>()))
        {
            return !(*this == other);
        }

        bool operator==(const std::vector<T>& other) const noexcept(noexcept(std::declval<T>() != std::declval<T>()))
        {
            if (this->size() != other.size()) return false;

            for (size_t col = 0; col < other.size(); col++)
            {
                if (mat_(row_, col) != other[col]) return false;
            }

            return true;
        }

        bool operator!=(const std::vector<T>& other) const noexcept(noexcept(std::declval<T>() != std::declval<T>()))
        {
            return !(*this == other);
        }
    };

    template<typename T, typename A>
    inline void swap(typename Matrix<T, A>::Row& lhs, typename Matrix<T, A>::Row& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    inline void swap(typename Matrix<T, A>::Row& lhs, std::vector<T, A>& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    inline void swap(std::vector<T, A>& lhs, typename Matrix<T, A>::Row& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        rhs.swap(lhs);
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_MATRIX_HPP