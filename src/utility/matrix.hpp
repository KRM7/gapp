/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATRIX_HPP
#define GA_UTILITY_MATRIX_HPP

#include "small_vector.hpp"
#include "iterators.hpp"
#include "functional.hpp"
#include "type_traits.hpp"
#include "utility.hpp"
#include <algorithm>
#include <vector>
#include <span>
#include <memory>
#include <initializer_list>
#include <type_traits>
#include <utility>
#include <cstddef>

namespace gapp::detail
{
    template<typename T, typename A>
    class MatrixRowBase;

    template<typename T, typename A>
    class MatrixRowRef;

    template<typename T, typename A>
    class ConstMatrixRowRef;


    template<typename T, typename A = std::allocator<T>>
    class Matrix : public iterator_interface<Matrix<T, A>>
    {
    public:
        using RowRef      = MatrixRowRef<T, A>;
        using ConstRowRef = ConstMatrixRowRef<T, A>;

        friend class MatrixRowBase<RowRef, Matrix>;
        friend class MatrixRowBase<ConstRowRef, const Matrix>;

        using value_type        = std::vector<T, A>;
        using allocator_type    = A;
        using storage_type      = std::vector<T, A>;

        using size_type         = typename storage_type::size_type;
        using difference_type   = typename storage_type::difference_type;

        using reference         = RowRef;
        using const_reference   = ConstRowRef;

        using element_type              = T;
        using element_reference         = typename storage_type::reference;
        using const_element_reference   = typename storage_type::const_reference;

        /* Iterators */

        class ConstRowIterator : public stable_iterator_base<ConstRowIterator, const Matrix, value_type, ConstRowRef, detail::proxy_ptr<ConstRowRef>, difference_type>
        {
        public:
            using my_base_ = stable_iterator_base<ConstRowIterator, const Matrix, value_type, ConstRowRef, detail::proxy_ptr<ConstRowRef>, difference_type>;
            using my_base_::my_base_;
            using typename my_base_::iterator_category;
        };

        class RowIterator : public stable_iterator_base<RowIterator, Matrix, value_type, RowRef, detail::proxy_ptr<RowRef>, difference_type>
        {
        public:
            using my_base_ = stable_iterator_base<RowIterator, Matrix, value_type, RowRef, detail::proxy_ptr<RowRef>, difference_type>;
            using my_base_::my_base_;
            using typename my_base_::iterator_category;

            constexpr /* implicit */ operator ConstRowIterator() const noexcept { return { this->data_, this->idx_ }; }
        };

        using iterator               = RowIterator;
        using const_iterator         = ConstRowIterator;
        using reverse_iterator       = std::reverse_iterator<RowIterator>;
        using const_reverse_iterator = std::reverse_iterator<ConstRowIterator>;

        constexpr iterator begin() noexcept { return iterator(*this, 0); }
        constexpr iterator end() noexcept   { return iterator(*this, nrows()); }

        constexpr const_iterator begin() const noexcept { return const_iterator(*this, 0); }
        constexpr const_iterator end() const noexcept   { return const_iterator(*this, nrows()); }

        /* Special member functions */

        Matrix()                         = default;
        Matrix(const Matrix&)            = default;
        Matrix(Matrix&&)                 = default;
        Matrix& operator=(const Matrix&) = default;
        Matrix& operator=(Matrix&&)      = default;
        ~Matrix()                        = default;

        constexpr explicit Matrix(const A& alloc) :
            data_(alloc), ncols_(0)
        {}

        constexpr Matrix(size_type nrows, size_type ncols) :
            data_(nrows * ncols), ncols_(ncols)
        {}

        constexpr Matrix(size_type nrows, size_type ncols, const T& init, const A& alloc = A()) :
            data_(nrows * ncols, init, alloc), ncols_(ncols)
        {}

        constexpr Matrix(std::initializer_list<std::initializer_list<T>> mat);

        constexpr Matrix(const_iterator first, const_iterator last);

        /* Member access */

        constexpr RowRef operator[](size_type row) noexcept
        {
            GAPP_ASSERT(row < nrows(), "Row index out of bounds.");

            return RowRef(*this, row);
        }

        constexpr ConstRowRef operator[](size_type row) const noexcept
        {
            GAPP_ASSERT(row < nrows(), "Row index out of bounds.");

            return ConstRowRef(*this, row);
        }

        constexpr storage_type column(size_type col_idx) const
        {
            GAPP_ASSERT(col_idx < ncols(), "Col index out of bounds.");

            storage_type col;
            col.reserve(nrows());

            for (size_type row = 0; row < nrows(); row++)
            {
                col.push_back((*this)(row, col_idx));
            }

            return col;
        }

        constexpr element_reference operator()(size_type row, size_type col) noexcept
        {
            GAPP_ASSERT(row < nrows(), "Row index out of bounds.");
            GAPP_ASSERT(col < ncols(), "Col index out of bounds.");

            return data_[row * ncols() + col];
        }

        constexpr const_element_reference operator()(size_type row, size_type col) const noexcept
        {
            GAPP_ASSERT(row < nrows(), "Row index out of bounds.");
            GAPP_ASSERT(col < ncols(), "Col index out of bounds.");

            return data_[row * ncols() + col];
        }

        constexpr RowRef front() noexcept { return *begin(); }
        constexpr ConstRowRef front() const noexcept { return *begin(); }

        constexpr RowRef back() noexcept { return *std::prev(end()); }
        constexpr ConstRowRef back() const noexcept { return *std::prev(end()); }

        constexpr auto data() noexcept { return data_.data(); }
        constexpr auto data() const noexcept { return data_.data(); }

        /* Modifiers (for rows) */

        constexpr void append_row(std::span<const T> row);
        constexpr void pop_back() noexcept { resize(nrows() - 1, ncols()); } // NOLINT(*exception-escape)

        /* Size / capacity */

        constexpr size_type size() const noexcept  { return nrows(); }
        constexpr size_type nrows() const noexcept { return empty() ? 0 : (data_.size() / ncols_); }
        constexpr size_type ncols() const noexcept { return ncols_; }
        constexpr bool empty() const noexcept { return data_.empty(); }

        constexpr void reserve(size_type nrows, size_type ncols) { data_.reserve(nrows * ncols); }

        constexpr void resize(size_type nrows, size_type ncols, const T& val = {})
        {
            data_.resize(nrows * ncols, val);
            ncols_ = ncols;
        }

        constexpr void clear() noexcept
        {
            data_.clear();
            ncols_ = 0;
        }

        /* Other */

        friend constexpr bool operator==(const Matrix& lhs, const Matrix& rhs)
        {
            return (lhs.empty() && rhs.empty()) ||
                   (lhs.nrows() == rhs.nrows() && lhs.ncols() == rhs.ncols() &&
                    lhs.data_ == rhs.data_);
        }

        constexpr void swap(Matrix& other) noexcept
        {
            std::swap(data_, other.data_);
            std::swap(ncols_, other.ncols_);
        }

    private:
        storage_type data_;
        size_type ncols_ = 0;
    };

    template<typename T, typename A>
    constexpr void swap(Matrix<T, A>& lhs, Matrix<T, A>& rhs) noexcept
    {
        lhs.swap(rhs);
    }


    /* MATRIX ROW PROXY CLASSES */

    template<typename Derived, typename MatrixType>
    class MatrixRowBase : public container_interface<Derived>
    {
    public:
        using value_type      = typename MatrixType::storage_type::value_type;
        using size_type       = typename MatrixType::storage_type::size_type;
        using difference_type = typename MatrixType::storage_type::difference_type;

        using reference = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::storage_type::const_reference,
                                                                          typename MatrixType::storage_type::reference>;
        using pointer   = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::storage_type::const_pointer,
                                                                          typename MatrixType::storage_type::pointer>;

        using const_reference = typename MatrixType::storage_type::const_reference;
        using const_pointer   = typename MatrixType::storage_type::const_pointer;

        using iterator = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::storage_type::const_iterator,
                                                                         typename MatrixType::storage_type::iterator>;
        using reverse_iterator = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::storage_type::const_reverse_iterator,
                                                                                 typename MatrixType::storage_type::reverse_iterator>;

        using const_iterator         = typename MatrixType::storage_type::const_iterator;
        using const_reverse_iterator = typename MatrixType::storage_type::const_reverse_iterator;

        constexpr MatrixRowBase(MatrixType& mat, size_type row) noexcept :
            mat_(&mat), row_(row)
        {}

        constexpr iterator begin() const noexcept
        {
            return mat_->data_.begin() + row_ * mat_->ncols();
        }

        constexpr iterator end() const noexcept
        {
            return begin() + mat_->ncols();
        }

        constexpr size_type ncols() const noexcept { return mat_->ncols(); }

        template<typename alloc_type>
        operator std::vector<value_type, alloc_type>() const { return std::vector<value_type, alloc_type>(begin(), end()); }

        operator small_vector<value_type>() const { return small_vector<value_type>(begin(), end()); }

        operator std::span<copy_const_t<MatrixType, value_type>>() const { return { begin(), end() }; }

        constexpr friend bool operator==(const Derived& lhs, const Derived& rhs)
        {
            if (lhs.size() != rhs.size()) return false;
            if (lhs.mat_ == rhs.mat_ && lhs.row_ == rhs.row_) return true;

            for (size_t col = 0; col < lhs.size(); col++)
            {
                if (lhs[col] != rhs[col]) return false;
            }

            return true;
        }

        template<typename alloc_type>
        constexpr friend bool operator==(const Derived& lhs, const std::vector<value_type, alloc_type>& rhs)
        {
            if (lhs.size() != rhs.size()) return false;

            for (size_t col = 0; col < lhs.size(); col++)
            {
                if (lhs[col] != rhs[col]) return false;
            }

            return true;
        }

    protected:
        MatrixType* mat_;
        size_type row_;
    };

    template<typename T, typename A>
    class MatrixRowRef : public MatrixRowBase<MatrixRowRef<T, A>, Matrix<T, A>>
    {
    public:
        friend class ConstMatrixRowRef<T, A>;

        using my_base_ = MatrixRowBase<MatrixRowRef<T, A>, Matrix<T, A>>;
        using my_base_::my_base_;

        using typename my_base_::const_iterator;

        MatrixRowRef(const MatrixRowRef&) = default;
        MatrixRowRef(MatrixRowRef&&)      = default;

        // NOLINTBEGIN(*unconventional-assign-operator, *assignment-signature)

        constexpr const MatrixRowRef& operator=(MatrixRowRef<T, A> rhs) const
        {
            return *this = std::span{ rhs.begin(), rhs.end() };
        }

        constexpr const MatrixRowRef& operator=(ConstMatrixRowRef<T, A> rhs) const
        {
            return *this = std::span{ rhs.begin(), rhs.end() };
        }

        constexpr const MatrixRowRef& operator=(std::initializer_list<T> rhs) const
        {
            return *this = std::span{ rhs.begin(), rhs.end() };
        }

        constexpr const MatrixRowRef& operator=(std::span<const T> rhs) const
        {
            GAPP_ASSERT(rhs.size() == this->size(), "Can't assign row with different length.");

            for (size_t col = 0; col < rhs.size(); col++)
            {
                (*this->mat_)(this->row_, col) = rhs[col];
            }

            return *this;
        }

        // NOLINTEND(*unconventional-assign-operator, *assignment-signature)

        constexpr void swap(std::span<T> other) const noexcept(std::is_nothrow_swappable_v<T>)
        {
            GAPP_ASSERT(other.size() == this->size(), "Rows must be the same size to swap them.");

            for (size_t col = 0; col < other.size(); col++)
            {
                using std::swap;
                swap((*this->mat_)(this->row_, col), other[col]);
            }
        }
    };

    template<typename T, typename A>
    constexpr void swap(MatrixRowRef<T, A> lhs, MatrixRowRef<T, A> rhs)
    noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    constexpr void swap(MatrixRowRef<T, A> lhs, std::span<std::type_identity_t<T>> rhs)
    noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    constexpr void swap(std::span<std::type_identity_t<T>> lhs, MatrixRowRef<T, A> rhs)
    noexcept(std::is_nothrow_swappable_v<T>)
    {
        rhs.swap(lhs);
    }


    template<typename T, typename A>
    class ConstMatrixRowRef : public MatrixRowBase<ConstMatrixRowRef<T, A>, const Matrix<T, A>>
    {
    public:
        friend class MatrixRowRef<T, A>;

        using my_base_ = MatrixRowBase<ConstMatrixRowRef<T, A>, const Matrix<T, A>>;
        using my_base_::my_base_;

        constexpr /* implicit */ ConstMatrixRowRef(MatrixRowRef<T, A> row) noexcept :
            my_base_(*row.mat_, row.row_)
        {}

        ConstMatrixRowRef(const ConstMatrixRowRef&) = default;
        ConstMatrixRowRef(ConstMatrixRowRef&&)      = default;

        ConstMatrixRowRef& operator=(const ConstMatrixRowRef&) = delete;
        ConstMatrixRowRef& operator=(ConstMatrixRowRef&&)      = delete;
    };


    /* MATRIX IMPLEMENTATION */

    template<typename T, typename A>
    constexpr Matrix<T, A>::Matrix(std::initializer_list<std::initializer_list<T>> mat)
    {
        if (mat.size() == 0) return;

        ncols_ = mat.begin()->size();

        GAPP_ASSERT(std::all_of(mat.begin(), mat.end(), detail::is_size(ncols_)), "Unequal row sizes in the input matrix.");

        data_.reserve(mat.size() * ncols_);
        for (auto& row : mat)
        {
            for (auto& entry : row)
            {
                data_.push_back(std::move(entry));
            }
        }
    }

    template<typename T, typename A>
    constexpr Matrix<T, A>::Matrix(const_iterator first, const_iterator last)
    {
        if (std::distance(first, last) <= 0) return;

        ncols_ = first->size();
        data_ = storage_type(first->begin(), first->begin() + std::distance(first, last) * ncols_);
    }

    template<typename T, typename A>
    constexpr void Matrix<T, A>::append_row(std::span<const T> row)
    {
        GAPP_ASSERT(row.size() == ncols() || nrows() == 0, "Can't insert row with different column count.");

        data_.insert(data_.end(), row.begin(), row.end());
        ncols_ = row.size();
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_MATRIX_HPP