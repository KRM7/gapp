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

        bool operator==(const Matrix&) const = default;
        bool operator!=(const Matrix&) const = default;

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

        /* Row iterators */

        iterator row_begin(size_type row) noexcept              { assert(row < nrows_); return data_.begin() + row * ncols_; }
        const_iterator row_begin(size_type row) const noexcept  { assert(row < nrows_); return data_.begin() + row * ncols_; }
        const_iterator row_cbegin(size_type row) const noexcept { assert(row < nrows_); return data_.cbegin() + row * ncols_; }
        iterator row_end(size_type row) noexcept                { assert(row < nrows_); return data_.begin() + (row + 1) * ncols_; }
        const_iterator row_end(size_type row) const noexcept    { assert(row < nrows_); return data_.begin() + (row + 1) * ncols_; }
        const_iterator row_cend(size_type row) const noexcept   { assert(row < nrows_); return data_.cbegin() + (row + 1) * ncols_; }

        reverse_iterator row_rbegin(size_type row) noexcept              { assert(row < nrows_); return data_.rbegin() + row * ncols_; }
        const_reverse_iterator row_rbegin(size_type row) const noexcept  { assert(row < nrows_); return data_.rbegin() + row * ncols_; }
        const_reverse_iterator row_crbegin(size_type row) const noexcept { assert(row < nrows_); return data_.crbegin() + row * ncols_; }
        reverse_iterator row_rend(size_type row) noexcept                { assert(row < nrows_); return data_.rbegin() + (row + 1) * ncols_; }
        const_reverse_iterator row_rend(size_type row) const noexcept    { assert(row < nrows_); return data_.rbegin() + (row + 1) * ncols_; }
        const_reverse_iterator row_crend(size_type row) const noexcept   { assert(row < nrows_); return data_.crbegin() + (row + 1) * ncols_; }

        /* Column iterators */

        class col_iterator
        {
        public:
            using iterator_category = std::iterator_traits<typename std::vector<T>::iterator>::iterator_category;
            using value_type        = std::iterator_traits<typename std::vector<T>::iterator>::value_type;
            using difference_type   = std::iterator_traits<typename std::vector<T>::iterator>::difference_type;
            using pointer           = std::iterator_traits<typename std::vector<T>::iterator>::pointer;
            using reference         = std::iterator_traits<typename std::vector<T>::iterator>::reference;

            col_iterator(std::vector<T>::iterator iter, difference_type stride) :
                iter_(iter), stride_(stride) {}

            bool operator==(const col_iterator& rhs) const noexcept { return iter_ == rhs.iter_; }
            bool operator!=(const col_iterator& rhs) const noexcept { return iter_ != rhs.iter_; }
            bool operator<(const col_iterator& rhs) const noexcept  { return iter_ < rhs.iter_; }
            bool operator<=(const col_iterator& rhs) const noexcept { return iter_ <= rhs.iter_; }
            bool operator>(const col_iterator& rhs) const noexcept  { return iter_ > rhs.iter_; }
            bool operator>=(const col_iterator& rhs) const noexcept { return iter_ >= rhs.iter_; }

            col_iterator& operator++() noexcept { iter_ += stride_; return *this; }
            col_iterator& operator--() noexcept { iter_ -= stride_; return *this; }
            col_iterator operator++(int) noexcept { auto temp = *this; iter_ += stride_; return temp; }
            col_iterator operator--(int) noexcept { auto temp = *this; iter_ -= stride_; return temp; }

            col_iterator& operator+=(difference_type n) noexcept { iter_ += (stride_ * n); return *this; }
            col_iterator& operator-=(difference_type n) noexcept { iter_ -= (stride_ * n); return *this; }
            col_iterator operator+(difference_type n) const { auto temp = col_iterator(*this); return temp += n; }
            col_iterator operator-(difference_type n) const { auto temp = col_iterator(*this); return temp -= n; }
            friend col_iterator operator+(difference_type, const col_iterator&) noexcept;

            difference_type operator-(const col_iterator& rhs) const noexcept { return (iter_ - rhs.iter_) % stride_; }

            reference operator*() const noexcept { return *iter_; }
            pointer operator->() const noexcept { return iter_.operator->(); }
            reference operator[](difference_type n) const noexcept { return iter_[stride_ * n]; }

        private:
            const difference_type stride_;
            std::vector<T>::iterator iter_;
        };


        class const_col_iterator
        {
        public:
            using iterator_category = std::iterator_traits<typename std::vector<T>::const_iterator>::iterator_category;
            using value_type        = std::iterator_traits<typename std::vector<T>::const_iterator>::value_type;
            using difference_type   = std::iterator_traits<typename std::vector<T>::const_iterator>::difference_type;
            using pointer           = std::iterator_traits<typename std::vector<T>::const_iterator>::pointer;
            using reference         = std::iterator_traits<typename std::vector<T>::const_iterator>::reference;

            const_col_iterator(std::vector<T>::const_iterator iter, difference_type stride) :
                iter_(iter), stride_(stride) {}

            const_col_iterator(std::vector<T>::iterator iter, difference_type stride) :
                iter_(iter), stride_(stride) {}

            bool operator==(const const_col_iterator& rhs) const noexcept { return iter_ == rhs.iter_; }
            bool operator!=(const const_col_iterator& rhs) const noexcept { return iter_ != rhs.iter_; }
            bool operator<(const const_col_iterator& rhs) const noexcept  { return iter_ < rhs.iter_; }
            bool operator<=(const const_col_iterator& rhs) const noexcept { return iter_ <= rhs.iter_; }
            bool operator>(const const_col_iterator& rhs) const noexcept  { return iter_ > rhs.iter_; }
            bool operator>=(const const_col_iterator& rhs) const noexcept { return iter_ >= rhs.iter_; }

            const_col_iterator& operator++() noexcept { iter_ += stride_; return *this; }
            const_col_iterator& operator--() noexcept { iter_ -= stride_; return *this; }
            const_col_iterator operator++(int) noexcept { auto temp = *this; iter_ += stride_; return temp; }
            const_col_iterator operator--(int) noexcept { auto temp = *this; iter_ -= stride_; return temp; }

            const_col_iterator& operator+=(difference_type n) noexcept { iter_ += (stride_ * n); return *this; }
            const_col_iterator& operator-=(difference_type n) noexcept { iter_ -= (stride_ * n); return *this; }
            const_col_iterator operator+(difference_type n) const { auto temp = const_col_iterator(*this); return temp += n; }
            const_col_iterator operator-(difference_type n) const { auto temp = const_col_iterator(*this); return temp -= n; }
            friend const_col_iterator operator+(difference_type, const const_col_iterator&) noexcept;

            difference_type operator-(const const_col_iterator& rhs) const noexcept { return (iter_ - rhs.iter_) % stride_; }

            reference operator*() const noexcept { return *iter_; }
            pointer operator->() const noexcept { return iter_.operator->(); }
            reference operator[](difference_type n) const noexcept { return iter_[stride_ * n]; }

        private:
            const difference_type stride_;
            std::vector<T>::const_iterator iter_;
        };

        using reverse_col_iterator       = std::reverse_iterator<col_iterator>;
        using const_reverse_col_iterator = std::reverse_iterator<const_col_iterator>;


        col_iterator col_begin(size_type col) noexcept              { assert(col < ncols_); return { data_.begin() + col, ncols_ }; }
        const_col_iterator col_begin(size_type col) const noexcept  { assert(col < ncols_); return { data_.begin() + col, ncols_ }; }
        const_col_iterator col_cbegin(size_type col) const noexcept { assert(col < ncols_); return { data_.cbegin() + col, ncols_ }; }
        col_iterator col_end(size_type col) noexcept                { assert(col < ncols_); return { data_.end() + col, ncols_}; }
        const_col_iterator col_end(size_type col) const noexcept    { assert(col < ncols_); return { data_.end() + col, ncols_ }; }
        const_col_iterator col_cend(size_type col) const noexcept   { assert(col < ncols_); return { data_.cend() + col, ncols_ }; }

        reverse_col_iterator col_rbegin(size_type col) noexcept              { assert(col < ncols_); return std::make_reverse_iterator(col_iterator(data_.end() + col, ncols_)); }
        const_reverse_col_iterator col_rbegin(size_type col) const noexcept  { assert(col < ncols_); return std::make_reverse_iterator(const_col_iterator(data_.end() + col, ncols_)); }
        const_reverse_col_iterator col_crbegin(size_type col) const noexcept { assert(col < ncols_); return std::make_reverse_iterator(const_col_iterator(data_.cend() + col, ncols_)); }
        reverse_col_iterator col_rend(size_type col) noexcept                { assert(col < ncols_); return std::make_reverse_iterator(col_iterator(data_.begin() + col, ncols_)); }
        const_reverse_col_iterator col_rend(size_type col) const noexcept    { assert(col < ncols_); return std::make_reverse_iterator(const_col_iterator(data_.begin() + col, ncols_)); }
        const_reverse_col_iterator col_crend(size_type col) const noexcept   { assert(col < ncols_); return std::make_reverse_iterator(const_col_iterator(data_.cbegin() + col, ncols_)); }


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

        template<typename F>
        void for_each(F&& f)
        {
            for (size_t row = 0; row < nrows(); row++)
                for (size_t col = 0; col < ncols(); col++)
                    std::invoke(f, row, col, operator()(row, col));
        }

    private:
        size_type nrows_;
        size_type ncols_;
        std::vector<T> data_;
    };

    template<typename T>
    Matrix<T>::col_iterator operator+(typename Matrix<T>::col_iterator::difference_type n,
                                      const typename Matrix<T>::col_iterator& it) noexcept
    {
        return Matrix<T>::col_iterator(it) + n;
    }

    template<typename T>
    Matrix<T>::const_col_iterator operator+(typename Matrix<T>::const_col_iterator::difference_type n,
                                            const typename Matrix<T>::const_col_iterator& it) noexcept
    {
        return Matrix<T>::const_col_iterator(it) + n;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_MATRIX_HPP