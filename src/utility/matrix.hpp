/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATRIX_HPP
#define GA_UTILITY_MATRIX_HPP

#include <vector>
#include <iterator>
#include <type_traits>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<typename T>
    class Matrix
    {
    public:
        using value_type        = typename std::vector<T>::value_type;
        using allocator_type    = typename std::vector<T>::allocator_type;

        using size_type         = typename std::vector<T>::size_type;
        using difference_type   = typename std::vector<T>::difference_type;

        using reference         = typename std::vector<T>::reference;
        using pointer           = typename std::vector<T>::pointer;
        using const_reference   = typename std::vector<T>::const_reference;
        using const_pointer     = typename std::vector<T>::const_pointer;

        using iterator               = typename std::vector<T>::iterator;
        using reverse_iterator       = typename std::vector<T>::reverse_iterator;
        using const_iterator         = typename std::vector<T>::const_iterator;
        using const_reverse_iterator = typename std::vector<T>::const_reverse_iterator;

        /* Special member functions */

        Matrix()                         = default;
        Matrix(const Matrix&)            = default;
        Matrix(Matrix&&)                 = default;
        Matrix& operator=(const Matrix&) = default;
        Matrix& operator=(Matrix&&)      = default;
        ~Matrix()                        = default;

        explicit Matrix(size_t nrows) :
            nrows_(nrows), ncols_(nrows), data_(nrows * nrows)
        {}

        Matrix(size_t nrows, const T& init) :
            nrows_(nrows), ncols_(nrows), data_(nrows * nrows, init)
        {}

        Matrix(size_t nrows, size_t ncols) :
            nrows_(nrows), ncols_(ncols), data_(nrows * ncols)
        {}

        Matrix(size_t nrows, size_t ncols, const T& init) :
            nrows_(nrows), ncols_(ncols), data_(nrows * ncols, init)
        {}

        /* Operators */

        friend bool operator==(const Matrix&, const Matrix&) = default;

        /* Member access */

        reference operator()(size_type row, size_type col) noexcept { return data_[row * ncols_ + col]; }
        const_reference operator()(size_type row, size_type col) const noexcept { return data_[row * ncols_ + col]; }

        reference at(size_type row, size_type col) { return data_.at(row * ncols_ + col); }
        const_reference at(size_type row, size_type col) const { return data_.at(row * ncols_ + col); }

        pointer data() noexcept { return data_.data(); }
        const_pointer data() const noexcept { return data_.data(); }

        /* Size / capacity */

        size_type nrows() const noexcept { return nrows_; }
        size_type ncols() const noexcept { return ncols_; }
        size_type size() const noexcept { return data_.size(); }
        size_type max_size() const noexcept { return data_.max_size(); }
        size_type capacity() const noexcept { return data_.capacity(); }
        size_type empty() const noexcept { return data_.empty(); }

        void reserve(size_type new_capacity) { data_.reserve(new_capacity); }
        void reserve(size_type nrows, size_type ncols) { data_.reserve(nrows * ncols); }
        void resize(size_type nrows, size_type ncols) { data_.resize(nrows * ncols); nrows_ = nrows; ncols_ = ncols; }
        void shrink_to_fit() { data_.shrink_to_fit(); }
        void clear() noexcept { data_.clear(); }

        /* Iterators */

        iterator begin() noexcept                       { return data_.begin(); }
        const_iterator begin() const noexcept           { return data_.begin(); }
        const_iterator cbegin() const noexcept          { return data_.cbegin(); }
        iterator end() noexcept                         { return data_.end(); }
        const_iterator end() const noexcept             { return data_.end(); }
        const_iterator cend() const noexcept            { return data_.cend(); }

        reverse_iterator rbegin() noexcept              { return data_.rbegin(); }
        const_reverse_iterator rbegin() const noexcept  { return data_.rbegin(); }
        const_reverse_iterator crbegin() const noexcept { return data_.crbegin(); }
        reverse_iterator rend() noexcept                { return data_.rend(); }
        const_reverse_iterator rend() const noexcept    { return data_.rend(); }
        const_reverse_iterator crend() const noexcept   { return data_.crend(); }

        /* Other */

        void swap(Matrix& other) noexcept(std::is_nothrow_swappable_v<T>)
        {
            std::swap(data_, other.data_);
            std::swap(nrows_, other.nrows_);
            std::swap(ncols_, other.ncols_);
        }

        void transpose() noexcept(std::is_nothrow_swappable_v<T>)
        {
            for (size_type row = 0; row < nrows(); row++)
                for (size_type col = row + 1; col < ncols(); col++)
                    std::swap(operator()(row, col), operator()(col, row));
        }

    private:
        size_type nrows_;
        size_type ncols_;
        std::vector<T> data_;
    };

    template<typename T>
    void swap(Matrix<T>& lhs, Matrix<T>& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_MATRIX_HPP