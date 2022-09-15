/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ITERATORS_HPP
#define GA_UTILITY_ITERATORS_HPP

#include <iterator>
#include <cstddef>

namespace genetic_algorithm::detail
{
    /*
    * The following should be implemented in Derived:
    *   begin(), end() with const overloads
    */
    template<typename Derived>
    class iterator_interface
    {
    public:
        auto cbegin() const { return derived.begin(); };
        auto cend() const { return derived.end(); };
    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };

    /*
    * The following should be implemented in Derived:
    *   begin(), end() with const overloads
    */
    template<typename Derived>
    class reverse_iterator_interface : public iterator_interface<Derived>
    {
    public:
        auto rbegin()        { return std::make_reverse_iterator(derived().end()); }
        auto rbegin() const  { return std::make_reverse_iterator(derived().end()); }
        auto crbegin() const { return std::make_reverse_iterator(derived().end()); }
        auto rend()          { return std::make_reverse_iterator(derived().begin()); }
        auto rend() const    { return std::make_reverse_iterator(derived().begin()); }
        auto crend() const   { return std::make_reverse_iterator(derived().begin()); }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   prefix operator++
    */
    template<typename Derived>
    class input_iterator_interface
    {
    public:
        Derived operator++(int)
        {
            Derived old_value(derived());
            ++derived();
            return old_value;
        }

        auto operator->() const
        {
            return &*derived();
        }

        friend bool operator!=(const Derived& lhs, const Derived& rhs)
        {
            return !(lhs == rhs);
        }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   prefix operator++
    */
    template<typename Derived>
    using forward_iterator_interface = input_iterator_interface<Derived>;


    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   prefix operator++
    *   prefix operator--
    */
    template<typename Derived>
    class bidirectional_iterator_interface : public forward_iterator_interface<Derived>
    {
    public:
        Derived operator--(int)
        {
            Derived old_value(derived());
            --derived();
            return old_value;
        }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };

    
    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   operator<
    *   prefix operator++
    *   prefix operator--
    *   operator +=(n)
    *   operator-(it, it)
    */
    template<typename Derived>
    class random_access_iterator_interface : public bidirectional_iterator_interface<Derived>
    {
    public:
        using difference_type = std::ptrdiff_t;

        friend bool operator>(const Derived& lhs, const Derived& rhs)  { return rhs < lhs; }
        friend bool operator<=(const Derived& lhs, const Derived& rhs) { return !(rhs < lhs); }
        friend bool operator>=(const Derived& lhs, const Derived& rhs) { return !(lhs < rhs); }

        friend Derived operator+(const Derived& iter, difference_type n) { return iter += n; }
        friend Derived operator+(difference_type n, const Derived& iter) { return iter += n; }

        Derived& operator-=(difference_type n) { return *this += -n; }
        friend Derived operator-(const Derived& iter, difference_type n) { return iter -= n; }

        auto operator[](difference_type n) { return *(*this + n); }
        const auto operator[](difference_type n) const { return *(*this + n); }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /* Random access iterator that can't be invalidated. */
    template<typename Container>
    class stable_iterator : public random_access_iterator_interface<stable_iterator<Container>>
    {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = typename Container::value_type;
        using reference = typename Container::reference;
        using typename random_access_iterator_interface<stable_iterator>::difference_type;

        stable_iterator() :
            data_(nullptr), idx_(0)
        {}

        stable_iterator(Container& container, size_t idx) :
            data_(&container), idx_(idx)
        {}

        reference operator*() const
        {
            assert(data_ && data_->size() > idx_);
            return data_[idx_];
        }

        friend bool operator==(stable_iterator lhs, stable_iterator rhs)
        {
            assert(lhs.data_ == rhs.data_);
            return (lhs.data_ == rhs.data_) && (lhs.idx_ == rhs.idx_);
        }

        friend bool operator<(stable_iterator lhs, stable_iterator rhs)
        {
            assert(lhs.data_ == rhs.data_);
            return lhs.idx_ < rhs.idx_;
        }

        stable_iterator& operator++()
        {
            assert(data_);
            ++idx_;
            return *this;
        }

        stable_iterator& operator--()
        {
            assert(data_ && idx_ != 0);
            --idx_;
            return *this;
        }

        stable_iterator& operator+=(difference_type n)
        {
            assert(data_); assert(n < 0 ? (idx_ >= -n) : true);
            idx_ += n;
            return *this;
        }

        friend difference_type operator-(stable_iterator lhs, stable_iterator rhs)
        {
            assert(lhs.data_ && lhs.data_ == rhs.data_);
            return difference_type(lhs.idx_) - difference_type(rhs.idx_);
        }

    private:
        Container* data_;
        size_t idx_;
    };


} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP