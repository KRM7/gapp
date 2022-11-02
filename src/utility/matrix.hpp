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

namespace genetic_algorithm::detail
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
        using value_type     = T;
        using allocator_type = A;
        using storage_type   = std::vector<T, A>;

        using size_type         = typename storage_type::size_type;
        using difference_type   = typename storage_type::difference_type;

        using reference         = typename storage_type::reference;
        using pointer           = typename storage_type::pointer;
        using const_reference   = typename storage_type::const_reference;
        using const_pointer     = typename storage_type::const_pointer;

        using RowRef      = MatrixRowRef<T, A>;
        using ConstRowRef = ConstMatrixRowRef<T, A>;

        friend class MatrixRowBase<RowRef, Matrix>;
        friend class MatrixRowBase<ConstRowRef, const Matrix>;

        /* Iterators */

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

        explicit Matrix(const A& alloc) :
            data_(alloc), nrows_(0), ncols_(0)
        {}

        Matrix(size_t nrows, size_t ncols, const T& init = T(), const A& alloc = A()) :
            data_(nrows * ncols, init, alloc), nrows_(nrows), ncols_(ncols)
        {}

        Matrix(std::initializer_list<std::initializer_list<T>> mat);

        /* Member access */

        RowRef operator[](size_t row) noexcept
        {
            GA_ASSERT(row < nrows_, "Row index out of bounds.");

            return RowRef(*this, row);
        }

        ConstRowRef operator[](size_t row) const noexcept
        {
            GA_ASSERT(row < nrows_, "Row index out of bounds.");

            return ConstRowRef(*this, row);
        }

        reference operator()(size_type row, size_type col) noexcept
        {
            GA_ASSERT(row < nrows_, "Row index out of bounds.");
            GA_ASSERT(col < ncols_, "Col index out of bounds.");

            return data_[row * ncols_ + col];
        }

        const_reference operator()(size_type row, size_type col) const noexcept
        {
            GA_ASSERT(row < nrows_, "Row index out of bounds.");
            GA_ASSERT(col < ncols_, "Col index out of bounds.");

            return data_[row * ncols_ + col];
        }

        pointer data() noexcept             { return data_.data(); }
        const_pointer data() const noexcept { return data_.data(); }

        /* Modifiers (for rows) */

        void append_row(const std::vector<T, A>& row);
        void append_row(std::vector<T, A>&& row);
        void append_row(ConstRowRef row);

        iterator erase(const_iterator row);
        iterator erase(const_iterator first, const_iterator last);

        /* Size / capacity */

        size_type size() const noexcept  { return nrows_; } /* For the bounds checking in stable_iterator */
        size_type nrows() const noexcept { return nrows_; }
        size_type ncols() const noexcept { return ncols_; }
        size_type empty() const noexcept { return data_.empty(); }

        void reserve(size_type nrows, size_type ncols) { data_.reserve(nrows * ncols); }
        void resize(size_type nrows, size_type ncols, const T& val = {}) { data_.resize(nrows * ncols, val); nrows_ = nrows; ncols_ = ncols; }

        /* Other */

        friend bool operator==(const Matrix& lhs, const Matrix& rhs)
        {
            return (lhs.empty() && rhs.empty()) ||
                   (lhs.nrows_ == rhs.nrows_ &&
                    lhs.ncols_ == rhs.ncols_ &&
                    lhs.data_ == rhs.data_);
        }

        friend bool operator!=(const Matrix& lhs, const Matrix& rhs)
        {
            return !(lhs == rhs);
        }

        void swap(Matrix& other) noexcept(std::is_nothrow_swappable_v<T>)
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
            /* implicit */ PtrHelper(const std::remove_cvref_t<U>& row) : row_(row) {}
            U* operator->() { return &row_; }
            U row_;
        };
    };

    template<typename T, typename A>
    inline void swap(Matrix<T, A>& lhs, Matrix<T, A>& rhs) noexcept(std::is_nothrow_swappable_v<T>)
    {
        lhs.swap(rhs);
    }


    /* MATRIX ROW PROXY IMPLEMENTATION */

    template<typename Derived, typename MatrixType>
    class MatrixRowBase : public reverse_iterator_interface<Derived>
    {
    public:
        using value_type      = typename MatrixType::storage_type::value_type;
        using size_type       = typename MatrixType::storage_type::size_type;
        using difference_type = typename MatrixType::storage_type::difference_type;

        using reference = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::const_reference, typename MatrixType::reference>;
        using pointer   = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::const_pointer, typename MatrixType::pointer>;

        using const_reference = typename MatrixType::const_reference;
        using const_pointer   = typename MatrixType::const_pointer;

        using iterator = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::storage_type::const_iterator,
                                                                         typename MatrixType::storage_type::iterator>;
        using reverse_iterator = std::conditional_t<std::is_const_v<MatrixType>, typename MatrixType::storage_type::const_reverse_iterator,
                                                                                 typename MatrixType::storage_type::reverse_iterator>;

        using const_iterator         = typename MatrixType::storage_type::const_iterator;
        using const_reverse_iterator = typename MatrixType::storage_type::const_reverse_iterator;

        MatrixRowBase(MatrixType& mat, size_type row) noexcept :
            mat_(&mat), row_(row)
        {}

        iterator begin() const noexcept
        {
            return mat_->data_.begin() + row_ * mat_->ncols();
        }

        iterator end() const noexcept
        {
            return begin() + mat_->ncols();
        }

        reference operator[](size_type col) const noexcept
        {
            GA_ASSERT(mat_->ncols() > col, "Col index out of bounds.");

            return (*mat_)(row_, col);
        }

        size_type size() const noexcept  { return mat_->ncols(); }
        size_type ncols() const noexcept { return mat_->ncols(); }

        /* Comparison operators */

        friend bool operator==(const Derived& lhs, const Derived& rhs) noexcept
        {
            if (lhs.size() != rhs.size()) return false;

            if (lhs.mat_ == rhs.mat_ && lhs.row_ == rhs.row_)
            {
                return true;
            }

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

        friend bool operator!=(const Derived& lhs, const Derived& rhs) noexcept                 { return !(lhs == rhs); }
        friend bool operator!=(const Derived& lhs, const std::vector<value_type>& rhs) noexcept { return !(lhs == rhs); }
        friend bool operator==(const std::vector<value_type>& lhs, const Derived& rhs) noexcept { return rhs == lhs; }
        friend bool operator!=(const std::vector<value_type>& lhs, const Derived& rhs) noexcept { return rhs != lhs; }

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
            return this->mat_->data_.cbegin() + this->row_ * this->mat_->ncols();
        }

        const_iterator cend() const noexcept
        {
            return cbegin() + this->mat_->ncols();
        }

        const MatrixRowRef& operator=(MatrixRowRef<T, A> rhs) const
        {
            return *this = ConstMatrixRowRef(rhs);
        }

        const MatrixRowRef& operator=(ConstMatrixRowRef<T, A> rhs) const
        {
            if (rhs.size() != this->size())
            {
                GA_THROW(std::invalid_argument, "Can't assign row with different length.");
            }

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

        const MatrixRowRef& operator=(const std::vector<T>& rhs) const
        {
            if (rhs.size() != this->size())
            {
                GA_THROW(std::invalid_argument, "Can't assign row with different length.");
            }

            for (size_t col = 0; col < rhs.size(); col++)
            {
                (*this->mat_)(this->row_, col) = rhs[col];
            }

            return *this;
        }

        const MatrixRowRef& operator=(std::vector<T>&& rhs) const
        {
            if (rhs.size() != this->size())
            {
                GA_THROW(std::invalid_argument, "Can't assign row with different length.");
            }

            for (size_t col = 0; col < rhs.size(); col++)
            {
                (*this->mat_)(this->row_, col) = std::move_if_noexcept(rhs[col]);
            }

            return *this;
        }

        void swap(MatrixRowRef other) const
        {
            if (other.size() != this->size())
            {
                GA_THROW(std::invalid_argument, "Rows must be the same size to swap them.");
            }

            for (size_t col = 0; col < other.size(); col++)
            {
                using std::swap;
                swap((*this->mat_)(this->row_, col), other[col]);
            }
        }

        void swap(std::vector<T>& other) const
        {
            if (other.size() != this->size())
            {
                GA_THROW(std::invalid_argument, "Rows must be the same size to swap them.");
            }

            for (size_t col = 0; col < other.size(); col++)
            {
                using std::swap;
                swap((*this->mat_)(this->row_, col), other[col]);
            }
        }
    };

    template<typename T, typename A>
    inline void swap(MatrixRowRef<T, A> lhs, MatrixRowRef<T, A> rhs)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    inline void swap(MatrixRowRef<T, A> lhs, std::vector<T, A>& rhs)
    {
        lhs.swap(rhs);
    }

    template<typename T, typename A>
    inline void swap(std::vector<T, A>& lhs, MatrixRowRef<T, A> rhs)
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

        /* implicit */ ConstMatrixRowRef(MatrixRowRef<T, A> row) noexcept :
            my_base_(*row.mat_, row.row_)
        {}

        ConstMatrixRowRef(const ConstMatrixRowRef&) = default;
        ConstMatrixRowRef(ConstMatrixRowRef&&)      = default;

        ConstMatrixRowRef& operator=(const ConstMatrixRowRef&) = delete;
        ConstMatrixRowRef& operator=(ConstMatrixRowRef&&)      = delete;
    };
    

    /* MATRIX ROW ITERATORS IMPLEMENTATIONS */

    template<typename T, typename A>
    class Matrix<T, A>::RowIterator :
        public stable_iterator_base<typename Matrix<T, A>::RowIterator, Matrix<T, A>,
                                    typename Matrix<T, A>::RowRef, typename Matrix<T, A>::RowRef, typename Matrix<T, A>::RowRef,
                                    typename Matrix<T, A>::difference_type>
    {
    public:
        using my_base_ = stable_iterator_base<RowIterator, Matrix, RowRef, RowRef, RowRef, difference_type>;
        using my_base_::my_base_;

        PtrHelper<RowRef> operator->() const { return **this; }

        friend class ConstRowIterator;
    };


    template<typename T, typename A>
    class Matrix<T, A>::ConstRowIterator :
        public stable_iterator_base<typename Matrix<T, A>::ConstRowIterator, const Matrix<T, A>,
                                    typename Matrix<T, A>::ConstRowRef, typename Matrix<T, A>::ConstRowRef, typename Matrix<T, A>::ConstRowRef,
                                    typename Matrix<T, A>::difference_type>
    {
    public:
        using my_base_ = stable_iterator_base<ConstRowIterator, const Matrix, ConstRowRef, ConstRowRef, ConstRowRef, difference_type>;
        using my_base_::my_base_;

        /* implicit */ ConstRowIterator(RowIterator it) noexcept :
            my_base_(*it.data_, it.idx_)
        {}

        PtrHelper<ConstRowRef> operator->() const { return **this; }
    };


    /* MATRIX IMPLEMENTATION */

    template<typename T, typename A>
    Matrix<T, A>::Matrix(std::initializer_list<std::initializer_list<T>> mat) :
        data_(), nrows_(mat.size()), ncols_(0)
    {
        if (mat.size() == 0) return;

        ncols_ = mat.begin()->size();

        if (!std::all_of(mat.begin(), mat.end(), is_size(ncols_)))
        {
            GA_THROW(std::invalid_argument, "Unequal row sizes in the input matrix.");
        }

        data_.reserve(nrows_ * ncols_);
        for (auto& row : mat)
        {
            for (auto& entry : row)
            {
                data_.push_back(std::move_if_noexcept(entry));
            }
        }
    }

    template<typename T, typename A>
    void Matrix<T, A>::append_row(const std::vector<T, A>& row)
    {
        if (row.size() != ncols_ && nrows_ != 0)
        {
            GA_THROW(std::invalid_argument, "Can't insert row with different column count.");
        }

        data_.insert(data_.end(), row.begin(), row.end());
        if (nrows_ == 0) ncols_ = row.size();
        nrows_++;
    }

    template<typename T, typename A>
    void Matrix<T, A>::append_row(std::vector<T, A>&& row)
    {
        if (row.size() != ncols_ && nrows_ != 0)
        {
            GA_THROW(std::invalid_argument, "Can't insert row with different column count.");
        }

        data_.insert(data_.end(), std::make_move_iterator(row.begin()), std::make_move_iterator(row.end()));
        if (nrows_ == 0) ncols_ = row.size();
        nrows_++;
    }

    template<typename T, typename A>
    void Matrix<T, A>::append_row(ConstRowRef row)
    {
        if (row.size() != ncols_ && nrows_ != 0)
        {
            GA_THROW(std::invalid_argument, "Can't insert row with different column count.");
        }

        data_.insert(data_.end(), row.begin(), row.end());
        if (nrows_ == 0) ncols_ = row.size();
        nrows_++;
    }

    template<typename T, typename A>
    auto Matrix<T, A>::erase(const_iterator row) -> iterator
    {
        const auto last_removed = data_.erase(row->begin(), row->end());
        --nrows_;

        return begin() + std::distance(data_.begin(), last_removed) / ncols_;
    }

    template<typename T, typename A>
    auto Matrix<T, A>::erase(const_iterator first, const_iterator last) -> iterator
    {
        const auto last_removed = data_.erase(first->begin(), last->begin());
        nrows_ -= std::distance(first, last);

        return begin() + std::distance(data_.begin(), last_removed) / ncols_;
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_MATRIX_HPP