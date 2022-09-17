/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATRIX_HPP
#define GA_UTILITY_MATRIX_HPP

#include "iterators.hpp"
#include <ostream>
#include <vector>
#include <memory>
#include <type_traits>
#include <utility>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<typename T, typename A>
    class RowReference;

    template<typename T, typename A = std::allocator<T>>
    class Matrix : public reverse_iterator_interface<Matrix<T, A>>
    {
    public:
        using value_type     = T;
        using allocator_type = A;
        using storage_type   = std::vector<T, A>;

        using size_type         = typename storage_type::size_type;
        using difference_type   = typename storage_type::difference_type;

        using reference         = typename storage_type::reference;
        using pointer           = typename storage_type::pointer;
        using const_reference   = typename storage_type::const_reference;
        using const_pointer     = typename storage_type::const_pointer;

        /* Iterators */

        friend class RowReference<T, A>; /* Not nested for free swap template argument deduction */
        using RowReference = RowReference<T, A>;

        class RowIterator;
        class ConstRowIterator;

        using iterator               = RowIterator;
        using const_iterator         = ConstRowIterator;
        using reverse_iterator       = std::reverse_iterator<RowIterator>;
        using const_reverse_iterator = std::reverse_iterator<ConstRowIterator>;

        iterator begin() noexcept { return iterator(*this, 0); }
        iterator end() noexcept   { return iterator(*this, this->size()); }

        const_iterator begin() const noexcept { return const_iterator(*this, 0); }
        const_iterator end() const noexcept   { return const_iterator(*this, this->size()); }

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

        Matrix(size_t nrows, size_t ncols, const T& init = T(), const A& alloc = A()) :
            nrows_(nrows), ncols_(ncols), data_(nrows * ncols, init, alloc)
        {}

        /* Operators */

        friend bool operator==(const Matrix&, const Matrix&) = default;

        RowReference operator[](size_t row) noexcept
        {
            return RowReference(*this, row);
        }

        const RowReference operator[](size_t row) const noexcept
        {
            return RowReference(*this, row);
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

        void push_back(const RowReference& row)
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

        size_type size() const noexcept { return nrows_; } /* For the bounds checking in stable_iterator */
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

        friend std::ostream& operator<<(std::ostream& out, const Matrix& mat)
        {
            for (const auto& row : mat)
            {
                out << "\n[ ";
                for (const auto& elem : row)
                {
                    out << elem << ", ";
                }
                out << "\b\b],";
            }
            out << "\b\n";

            return out;
        }

    private:
        storage_type data_;
        size_type nrows_;
        size_type ncols_;

        /* Helper struct for the iterators' operator-> */
        template<typename T>
        struct PtrHelper
        {
            /* implicit */ PtrHelper(const std::remove_cvref_t<T>& row) : row_(row) {}
            T* operator->() const { return &row_; }
            T row_;
        };
    };

    template<typename T, typename A>
    inline void swap(Matrix<T, A>& lhs, Matrix<T, A>& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }


    /* MATRIX ROW PROXY IMPLEMENTATION */

    template<typename T, typename A>
    class RowReference : public reverse_iterator_interface<RowReference<T, A>>
    {
    public:
        using value_type        = typename Matrix<T, A>::storage_type::value_type;
        using size_type         = typename Matrix<T, A>::storage_type::size_type;
        using difference_type   = typename Matrix<T, A>::storage_type::difference_type;
        using reference         = typename Matrix<T, A>::storage_type::reference;
        using const_reference   = typename Matrix<T, A>::storage_type::const_reference;


        RowReference(const Matrix<T, A>& mat, size_t row) noexcept :
            mat_(const_cast<Matrix<T, A>&>(mat)), row_(row) { assert(mat.nrows() > row); }

        RowReference(const RowReference&) = default;
        RowReference(RowReference&&)      = default;

        ~RowReference() = default;


        using iterator               = typename Matrix<T, A>::storage_type::iterator;
        using const_iterator         = typename Matrix<T, A>::storage_type::const_iterator;
        using reverse_iterator       = typename Matrix<T, A>::storage_type::reverse_iterator;
        using const_reverse_iterator = typename Matrix<T, A>::storage_type::const_reverse_iterator;

        iterator begin() noexcept              { return mat_.data_.begin() + row_ * mat_.ncols(); }
        iterator end() noexcept                { return mat_.data_.begin() + (row_ + 1) * mat_.ncols(); }

        const_iterator begin() const noexcept  { return mat_.data_.cbegin() + row_ * mat_.ncols(); }
        const_iterator end() const noexcept    { return mat_.data_.cbegin() + (row_ + 1) * mat_.ncols(); }

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

        reference at(size_t col) { return mat_.at(row_, col); }
        const_reference at(size_t col) const { return mat_.at(row_, col); }

        size_type size() const noexcept { return mat_.ncols(); }


        RowReference& operator=(RowReference rhs) noexcept(std::is_nothrow_assignable_v<T, T>)
        {
            assert(rhs.size() == this->size());

            for (size_t col = 0; col < rhs.size(); col++)
            {
                mat_(row_, col) = rhs[col];
            }

            return *this;
        }

        RowReference& operator=(const std::vector<T>& rhs) noexcept(std::is_nothrow_assignable_v<T, T>)
        {
            assert(rhs.size() == this->size());

            for (size_t col = 0; col < rhs.size(); col++)
            {
                mat_(row_, col) = rhs[col];
            }

            return *this;
        }

        RowReference& operator=(std::vector<T>&& rhs)
        {
            assert(rhs.size() == this->size());

            for (size_t col = 0; col < rhs.size(); col++)
            {
                mat_(row_, col) = std::move_if_noexcept(rhs[col]);
            }

            return *this;
        }


        void swap(RowReference other) noexcept(std::is_nothrow_swappable_v<T>)
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

        bool operator==(const RowReference& other) const noexcept(noexcept(std::declval<T>() != std::declval<T>()))
        {
            if ((&mat_ == &other.mat_) && (row_ == other.row_)) return true;

            if (this->size() != other.size()) return false;

            for (size_t col = 0; col < other.size(); col++)
            {
                if (mat_(row_, col) != other[col]) return false;
            }

            return true;
        }

        bool operator!=(const RowReference& other) const noexcept(noexcept(std::declval<T>() != std::declval<T>()))
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

    private:
        Matrix<T, A>& mat_;
        size_t row_;
    };

    template<typename T, typename A>
    inline void swap(RowReference<T, A> lhs, RowReference<T, A> rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    inline void swap(RowReference<T, A> lhs, std::vector<T, A>& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    inline void swap(std::vector<T, A>& lhs, RowReference<T, A> rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        rhs.swap(lhs);
    }
    

    /* MATRIX ROW ITERATORS IMPLEMENTATIONS */

    template<typename T, typename A>
    class Matrix<T, A>::RowIterator :
        public stable_iterator_base<typename Matrix<T, A>::RowIterator,
                                    typename Matrix<T, A>::RowReference, typename Matrix<T, A>::RowReference, typename Matrix<T, A>::RowReference,
                                    typename Matrix<T, A>::difference_type>
    {
    public:
        using _my_base = stable_iterator_base<RowIterator, RowReference, RowReference, RowReference, difference_type>;

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using typename _my_base::value_type;
        using typename _my_base::reference;
        using typename _my_base::pointer;

        RowIterator() noexcept :
            _my_base(), data_(nullptr), idx_(0)
        {}

        RowIterator(Matrix& mat, size_t row) noexcept :
            _my_base(), data_(&mat), idx_(row)
        {}

        PtrHelper<RowReference> operator->() const
        {
            return (*data_)[idx_];
        }

    private:
        Matrix* data_;
        size_t idx_;

        friend class _my_base;
    };

    template<typename T, typename A>
    class Matrix<T, A>::ConstRowIterator :
        public stable_iterator_base<typename Matrix<T, A>::ConstRowIterator,
                                    const typename Matrix<T, A>::RowReference, const typename Matrix<T, A>::RowReference, const typename Matrix<T, A>::RowReference,
                                    typename Matrix<T, A>::difference_type>
    {
    public:
        using _my_base = stable_iterator_base<ConstRowIterator, const RowReference, const RowReference, const RowReference, difference_type>;

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using typename _my_base::value_type;
        using typename _my_base::reference;
        using typename _my_base::pointer;

        ConstRowIterator() noexcept :
            _my_base(), data_(nullptr), idx_(0)
        {}

        ConstRowIterator(const Matrix& mat, size_t row) noexcept :
            _my_base(), data_(&mat), idx_(row)
        {}

        /* implicit */ ConstRowIterator(RowIterator it) noexcept :
            _my_base(), data_(it.data_), idx_(it.idx_)
        {}

        PtrHelper<const RowReference> operator->() const
        {
            return (*data_)[idx_];
        }

    private:
        const Matrix* data_;
        size_t idx_;

        friend class _my_base;
    };

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_MATRIX_HPP