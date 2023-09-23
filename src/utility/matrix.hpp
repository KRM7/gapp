/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_MATRIX_HPP
#define GA_UTILITY_MATRIX_HPP

#include "iterators.hpp"
#include "functional.hpp"
#include "utility.hpp"
#include <vector>
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
    class Matrix : public reverse_iterator_interface<Matrix<T, A>>
    {
    public:
        using RowRef      = MatrixRowRef<T, A>;
        using ConstRowRef = ConstMatrixRowRef<T, A>;

        friend class MatrixRowBase<RowRef, Matrix>;
        friend class MatrixRowBase<ConstRowRef, const Matrix>;

        class RowIterator;
        class ConstRowIterator;

        using value_type      = RowRef;
        using allocator_type  = A;
        using storage_type    = std::vector<T, A>;

        using size_type       = typename storage_type::size_type;
        using difference_type = typename storage_type::difference_type;

        using reference       = RowRef;
        using pointer         = RowIterator;
        using const_reference = ConstRowRef;
        using const_pointer   = ConstRowIterator;

        using element_type            = T;
        using element_reference       = typename storage_type::reference;
        using const_element_reference = typename storage_type::const_reference;
        using element_pointer         = typename storage_type::pointer;
        using const_element_pointer   = typename storage_type::const_pointer;

        /* Iterators */

        using iterator               = RowIterator;
        using const_iterator         = ConstRowIterator;
        using reverse_iterator       = std::reverse_iterator<RowIterator>;
        using const_reverse_iterator = std::reverse_iterator<ConstRowIterator>;

        constexpr iterator begin() noexcept { return iterator(*this, 0); }
        constexpr iterator end() noexcept   { return iterator(*this, this->size()); }

        constexpr const_iterator begin() const noexcept { return const_iterator(*this, 0); }
        constexpr const_iterator end() const noexcept   { return const_iterator(*this, this->size()); }

        /* Special member functions */

        Matrix()                         = default;
        Matrix(const Matrix&)            = default;
        Matrix(Matrix&&)                 = default;
        Matrix& operator=(const Matrix&) = default;
        Matrix& operator=(Matrix&&)      = default;
        ~Matrix()                        = default;

        constexpr explicit Matrix(const A& alloc) :
            data_(alloc), nrows_(0), ncols_(0)
        {}

        constexpr Matrix(size_t nrows, size_t ncols) :
            data_(nrows * ncols), nrows_(nrows), ncols_(ncols)
        {}

        constexpr Matrix(size_t nrows, size_t ncols, const T& init, const A& alloc = A()) :
            data_(nrows * ncols, init, alloc), nrows_(nrows), ncols_(ncols)
        {}

        constexpr Matrix(std::initializer_list<std::initializer_list<T>> mat);

        constexpr Matrix(const_iterator first, const_iterator last);

        /* Member access */

        constexpr RowRef operator[](size_t row) noexcept
        {
            GAPP_ASSERT(row < nrows_, "Row index out of bounds.");

            return RowRef(*this, row);
        }

        constexpr ConstRowRef operator[](size_t row) const noexcept
        {
            GAPP_ASSERT(row < nrows_, "Row index out of bounds.");

            return ConstRowRef(*this, row);
        }

        constexpr element_reference operator()(size_type row, size_type col) noexcept
        {
            GAPP_ASSERT(row < nrows_, "Row index out of bounds.");
            GAPP_ASSERT(col < ncols_, "Col index out of bounds.");

            return data_[row * ncols_ + col];
        }

        constexpr const_element_reference operator()(size_type row, size_type col) const noexcept
        {
            GAPP_ASSERT(row < nrows_, "Row index out of bounds.");
            GAPP_ASSERT(col < ncols_, "Col index out of bounds.");

            return data_[row * ncols_ + col];
        }

        constexpr RowRef front() noexcept { return operator[](0); }
        constexpr ConstRowRef front() const noexcept { return operator[](0); }

        constexpr RowRef back() noexcept { return operator[](nrows_ - 1); }
        constexpr ConstRowRef back() const noexcept { return operator[](nrows_ - 1); }

        constexpr element_pointer data() noexcept             { return data_.data(); }
        constexpr const_element_pointer data() const noexcept { return data_.data(); }

        /* Modifiers (for rows) */

        constexpr void append_row(const std::vector<T, A>& row);
        constexpr void append_row(std::vector<T, A>&& row);
        constexpr void append_row(ConstRowRef row);

        constexpr iterator erase(const_iterator row);
        constexpr iterator erase(const_iterator first, const_iterator last);

        /* Size / capacity */

        constexpr size_type size() const noexcept  { return nrows_; } /* For the bounds checking in stable_iterator */
        constexpr size_type nrows() const noexcept { return nrows_; }
        constexpr size_type ncols() const noexcept { return ncols_; }
        constexpr size_type empty() const noexcept { return data_.empty(); }

        constexpr void reserve(size_type nrows, size_type ncols) { data_.reserve(nrows * ncols); }

        constexpr void resize(size_type nrows, size_type ncols, const T& val = {})
        {
            data_.resize(nrows * ncols, val);
            nrows_ = nrows;
            ncols_ = ncols;
        }

        constexpr void clear() noexcept
        {
            data_.clear();
            nrows_ = 0;
        }

        /* Other */

        friend constexpr bool operator==(const Matrix& lhs, const Matrix& rhs)
        {
            return (lhs.empty() && rhs.empty()) ||
                   (lhs.nrows_ == rhs.nrows_ && lhs.ncols_ == rhs.ncols_ &&
                    lhs.data_ == rhs.data_);
        }

        constexpr void swap(Matrix& other) noexcept
        {
            std::swap(data_, other.data_);
            std::swap(nrows_, other.nrows_);
            std::swap(ncols_, other.ncols_);
        }

    private:
        storage_type data_;
        size_type nrows_ = 0;
        size_type ncols_ = 0;

        /* Helper struct for the iterators' operator-> */
        template<typename U>
        struct PtrHelper
        {
            /* implicit */ constexpr PtrHelper(const std::remove_cvref_t<U>& row) : row_(row) {}
            constexpr U* operator->() noexcept { return &row_; }
            U row_;
        };
    };

    template<typename T, typename A>
    constexpr void swap(Matrix<T, A>& lhs, Matrix<T, A>& rhs) noexcept
    {
        lhs.swap(rhs);
    }


    /* MATRIX ROW PROXY CLASSES */

    template<typename Derived, typename MatrixType>
    class MatrixRowBase : public reverse_iterator_interface<Derived>
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

        constexpr reference operator[](size_type col) const noexcept
        {
            GAPP_ASSERT(mat_->ncols() > col, "Col index out of bounds.");

            return (*mat_)(row_, col);
        }

        constexpr size_type size() const noexcept  { return mat_->ncols(); }
        constexpr size_type ncols() const noexcept { return mat_->ncols(); }

        constexpr bool empty() const noexcept { return size() == 0_sz; }

        explicit operator std::vector<value_type>() { return std::vector(begin(), end()); }

        /* Comparison operators */

        friend bool operator==(const Derived& lhs, const Derived& rhs) noexcept
        {
            if (lhs.size() != rhs.size()) return false;
            if (lhs.mat_ == rhs.mat_ && lhs.row_ == rhs.row_) return true;

            for (size_t col = 0; col < lhs.size(); col++)
            {
                if (lhs[col] != rhs[col]) return false;
            }

            return true;
        }

        friend bool operator==(const Derived& lhs, const std::vector<value_type>& rhs) noexcept
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

        const_iterator cbegin() const noexcept
        {
            return this->mat_->cbegin() + this->row_ * this->ncols();
        }

        constexpr const_iterator cend() const noexcept
        {
            return cbegin() + this->ncols();
        }

        constexpr const MatrixRowRef& operator=(const std::vector<T>& rhs) const;
        constexpr const MatrixRowRef& operator=(std::vector<T>&& rhs) const;
        constexpr const MatrixRowRef& operator=(ConstMatrixRowRef<T, A> rhs) const;
        constexpr const MatrixRowRef& operator=(MatrixRowRef<T, A> rhs) const { return *this = ConstMatrixRowRef(rhs); }

        constexpr void swap(MatrixRowRef other) const;
        constexpr void swap(std::vector<T>& other) const;
    };

    template<typename T, typename A>
    constexpr void swap(MatrixRowRef<T, A> lhs, MatrixRowRef<T, A> rhs)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    constexpr void swap(MatrixRowRef<T, A> lhs, std::vector<T, A>& rhs)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    constexpr void swap(std::vector<T, A>& lhs, MatrixRowRef<T, A> rhs)
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
    

    /* MATRIX ROW ITERATORS */

    template<typename T, typename A>
    class Matrix<T, A>::RowIterator :
        public stable_iterator_base<typename Matrix<T, A>::RowIterator, Matrix<T, A>,
                                    typename Matrix<T, A>::RowRef, typename Matrix<T, A>::RowRef, typename Matrix<T, A>::RowRef,
                                    typename Matrix<T, A>::difference_type>
    {
    public:
        using iterator_category = std::contiguous_iterator_tag;
        using my_base_ = stable_iterator_base<RowIterator, Matrix, RowRef, RowRef, RowRef, difference_type>;
        using my_base_::my_base_;

        constexpr PtrHelper<RowRef> operator->() const { return **this; }

        friend class ConstRowIterator;
    };


    template<typename T, typename A>
    class Matrix<T, A>::ConstRowIterator :
        public stable_iterator_base<typename Matrix<T, A>::ConstRowIterator, const Matrix<T, A>,
                                    typename Matrix<T, A>::ConstRowRef, typename Matrix<T, A>::ConstRowRef, typename Matrix<T, A>::ConstRowRef,
                                    typename Matrix<T, A>::difference_type>
    {
    public:
        using iterator_category = std::contiguous_iterator_tag;
        using my_base_ = stable_iterator_base<ConstRowIterator, const Matrix, ConstRowRef, ConstRowRef, ConstRowRef, difference_type>;
        using my_base_::my_base_;

        constexpr /* implicit */ ConstRowIterator(RowIterator it) noexcept
        {
            this->data_ = it.data_;
            this->idx_  = it.idx_;
        }

        constexpr PtrHelper<ConstRowRef> operator->() const { return **this; }
    };



    /* MATRIX IMPLEMENTATION */

    template<typename T, typename A>
    constexpr Matrix<T, A>::Matrix(std::initializer_list<std::initializer_list<T>> mat) :
        data_(), nrows_(mat.size()), ncols_(0)
    {
        if (mat.size() == 0) return;

        ncols_ = mat.begin()->size();

        GAPP_ASSERT(std::all_of(mat.begin(), mat.end(), detail::is_size(ncols_)), "Unequal row sizes in the input matrix.");

        data_.reserve(nrows_ * ncols_);
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

        nrows_ = std::distance(first, last);
        ncols_ = first->size();
        data_ = storage_type(first->begin(), first->begin() + nrows_ * ncols_);
    }

    template<typename T, typename A>
    constexpr void Matrix<T, A>::append_row(const std::vector<T, A>& row)
    {
        GAPP_ASSERT(row.size() == ncols_ || nrows_ == 0, "Can't insert row with different column count.");

        data_.insert(data_.end(), row.begin(), row.end());
        if (nrows_ == 0) ncols_ = row.size();
        nrows_++;
    }

    template<typename T, typename A>
    constexpr void Matrix<T, A>::append_row(std::vector<T, A>&& row)
    {
        GAPP_ASSERT(row.size() == ncols_ || nrows_ == 0, "Can't insert row with different column count.");

        data_.insert(data_.end(), std::make_move_iterator(row.begin()), std::make_move_iterator(row.end()));
        if (nrows_ == 0) ncols_ = row.size();
        nrows_++;
    }

    template<typename T, typename A>
    constexpr void Matrix<T, A>::append_row(ConstRowRef row)
    {
        GAPP_ASSERT(row.size() == ncols_ || nrows_ == 0, "Can't insert row with different column count.");

        data_.insert(data_.end(), row.begin(), row.end());
        if (nrows_ == 0) ncols_ = row.size();
        nrows_++;
    }

    template<typename T, typename A>
    constexpr auto Matrix<T, A>::erase(const_iterator row) -> iterator
    {
        const auto last_removed = data_.erase(row->begin(), row->end());
        --nrows_;

        return begin() + std::distance(data_.begin(), last_removed) / ncols_;
    }

    template<typename T, typename A>
    constexpr auto Matrix<T, A>::erase(const_iterator first, const_iterator last) -> iterator
    {
        const auto data_last = last != end() ? last->begin() : data_.end(); // can't dereference last iter

        const auto last_removed = data_.erase(first->begin(), data_last);
        nrows_ -= std::distance(first, last);

        return begin() + std::distance(data_.begin(), last_removed) / ncols_;
    }


    /* MATRIX ROW PROXY IMPLEMENTATIONS */

    template<typename T, typename A>
    constexpr const MatrixRowRef<T, A>& MatrixRowRef<T, A>::operator=(ConstMatrixRowRef<T, A> rhs) const
    {
        GAPP_ASSERT(rhs.size() == this->size(), "Can't assign row with different length.");

        if (rhs.mat_ == this->mat_ && rhs.row_ == this->row_)
        {
            return *this;
        }

        for (size_t col = 0; col < rhs.size(); col++)
        {
            (*this->mat_)(this->row_, col) = rhs[col];
        }

        return *this;
    }

    template<typename T, typename A>
    constexpr const MatrixRowRef<T, A>& MatrixRowRef<T, A>::operator=(const std::vector<T>& rhs) const
    {
        GAPP_ASSERT(rhs.size() == this->size(), "Can't assign row with different length.");

        for (size_t col = 0; col < rhs.size(); col++)
        {
            (*this->mat_)(this->row_, col) = rhs[col];
        }

        return *this;
    }

    template<typename T, typename A>
    constexpr const MatrixRowRef<T, A>& MatrixRowRef<T, A>::operator=(std::vector<T>&& rhs) const
    {
        GAPP_ASSERT(rhs.size() == this->size(), "Can't assign row with different length.");

        for (size_t col = 0; col < rhs.size(); col++)
        {
            (*this->mat_)(this->row_, col) = std::move(rhs[col]);
        }

        return *this;
    }

    template<typename T, typename A>
    constexpr void MatrixRowRef<T, A>::swap(MatrixRowRef other) const
    {
        GAPP_ASSERT(other.size() == this->size(), "Rows must be the same size to swap them.");

        for (size_t col = 0; col < other.size(); col++)
        {
            using std::swap;
            swap((*this->mat_)(this->row_, col), other[col]);
        }
    }

    template<typename T, typename A>
    constexpr void MatrixRowRef<T, A>::swap(std::vector<T>& other) const
    {
        GAPP_ASSERT(other.size() == this->size(), "Rows must be the same size to swap them.");

        for (size_t col = 0; col < other.size(); col++)
        {
            using std::swap;
            swap((*this->mat_)(this->row_, col), other[col]);
        }
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_MATRIX_HPP