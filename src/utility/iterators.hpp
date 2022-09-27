/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ITERATORS_HPP
#define GA_UTILITY_ITERATORS_HPP

#include "utility.hpp"
#include <iterator>
#include <type_traits>
#include <limits>
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
        auto cbegin() const { return derived().begin(); };
        auto cend() const { return derived().end(); };

    private:
        Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
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
        Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   increment() function (equivalent to prefix operator++)
    */
    template<typename Derived>
    class input_iterator_interface
    {
    public:
        using iterator_category = std::input_iterator_tag;

        Derived& operator++()
        {
            return derived().increment();
        }

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
        Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   increment() function (equivalent to prefix operator++)
    */
    template<typename Derived>
    class forward_iterator_interface : public input_iterator_interface<Derived>
    {
    public:
        using iterator_category = std::forward_iterator_tag;
    };


    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   increment() function (equivalent to prefix operator++)
    *   decrement() function (equivalent to prefix operator--)
    */
    template<typename Derived>
    class bidirectional_iterator_interface : public forward_iterator_interface<Derived>
    {
    public:
        using iterator_category = std::bidirectional_iterator_tag;

        Derived& operator--()
        {
            return derived().decrement();
        }

        Derived operator--(int)
        {
            Derived old_value(derived());
            --derived();
            return old_value;
        }

    private:
        Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };

    
    /*
    * The following should be implemented in Derived:
    *   operator*
    *   operator==
    *   operator<
    *   increment() function (equivalent to prefix operator++)
    *   decrement() function (equivalent to prefix operator--)
    *   operator +=(n)
    *   operator-(it, it)
    */
    template<typename Derived, typename Distance = std::ptrdiff_t>
    class random_access_iterator_interface : public bidirectional_iterator_interface<Derived>
    {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type = Distance;

        friend bool operator>(const Derived& lhs, const Derived& rhs)  { return rhs < lhs; }
        friend bool operator<=(const Derived& lhs, const Derived& rhs) { return !(rhs < lhs); }
        friend bool operator>=(const Derived& lhs, const Derived& rhs) { return !(lhs < rhs); }

        friend Derived operator+(Derived iter, difference_type n) { return iter += n; }
        friend Derived operator+(difference_type n, Derived iter) { return iter += n; }

        friend Derived& operator-=(Derived& lhs, difference_type n) { return lhs += -n; }
        friend Derived operator-(Derived iter, difference_type n) { return iter -= n; }

        decltype(auto) operator[](difference_type n) const { return *(derived() + n); }

    private:
        Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /* Iterators for random-access containers that aren't invalidated on reallocations and insertions. */

    template<typename Derived,
             typename Container,
             typename ValueType,
             typename Reference,
             typename Pointer,
             typename Distance = std::ptrdiff_t>
    class stable_iterator_base : public random_access_iterator_interface<Derived, Distance>
    {
    public:
        using _my_base = random_access_iterator_interface<Derived>;

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using value_type = ValueType;
        using reference  = Reference;
        using pointer    = Pointer;

        stable_iterator_base() noexcept :
            data_(nullptr), idx_(0)
        {}

        stable_iterator_base(Container& container, size_t idx) noexcept :
            data_(&container), idx_(idx)
        {
            GA_ASSERT(data_->size() >= idx_, "Iterator can't refer to element past the end of the range.");
        }

        reference operator*() const
        {
            GA_ASSERT(data_ != nullptr, "Can't dereference value initialized iterator.");
            GA_ASSERT(data_->size() > idx_, "Can't dereference past-the-end iterator.");

            return (*data_)[idx_];
        }

        bool operator==(const Derived& rhs) const
        {
            GA_ASSERT(data_ == rhs.data_, "Can't compare iterators of different ranges.");
            GA_ASSERT(data_ == nullptr || data_->size() >= idx_, "Can't compare invalid iterator.");
            GA_ASSERT(rhs.data_ == nullptr || rhs.data_->size() >= rhs.idx_, "Can't compare invalid iterator.");

            return idx_ == rhs.idx_;    /* Value-initialized iterators will have the same idx. */
        }

        bool operator<(const Derived& rhs) const
        {
            GA_ASSERT(data_ == rhs.data_, "Can't compare iterators of different ranges.");
            GA_ASSERT(data_ == nullptr || data_->size() >= idx_, "Can't compare invalid iterator.");
            GA_ASSERT(rhs.data_ == nullptr || rhs.data_->size() >= rhs.idx_, "Can't compare invalid iterator.");

            return idx_ < rhs.idx_;     /* Value-initialized iterators will have the same idx. */
        }

        Derived& increment() noexcept
        {
            GA_ASSERT(data_ != nullptr, "Can't increment value initialized iterator.");
            GA_ASSERT(idx_ != data_->size(), "Can't increment past-the-end iterator.");

            ++idx_;
            return static_cast<Derived&>(*this);
        }

        Derived& decrement() noexcept
        {
            GA_ASSERT(data_ != nullptr, "Can't decrement value initialized iterator.");
            GA_ASSERT(idx_ != 0, "Can't decremenet the begin iterator.");

            --idx_;
            return static_cast<Derived&>(*this);
        }

        Derived& operator+=(difference_type n)
        {
            GA_ASSERT(data_ != nullptr, "Can't offset value initialized iterator.");
            GA_ASSERT(n < 0 ? (difference_type(idx_) >= -n) : true, "Can't move iterator to before the start of the range.");
            GA_ASSERT(n > 0 ? (idx_ + n) <= data_->size() : true, "Can't move iterator past the end of the range.");

            idx_ += n;
            return static_cast<Derived&>(*this);
        }

        difference_type operator-(const Derived& rhs) const
        {
            GA_ASSERT(data_ && rhs.data_, "Can't get the distance of value initialized iterators.");
            GA_ASSERT(data_ == rhs.data_, "Can't get the distance of iterators of different ranges.");
            GA_ASSERT(data_->size() >= idx_, "Invalid lhs iterator.");
            GA_ASSERT(rhs.data_->size() >= rhs.idx_, "Invalid rhs iterator.");
            GA_ASSERT((idx_ >= rhs.idx_ ? idx_ - rhs.idx_ : rhs.idx_ - idx_) <= size_t(std::numeric_limits<difference_type>::max()),
                     "Can't represent the result of the operation as difference_type.");

            return idx_ >= rhs.idx_ ?
                difference_type(idx_ - rhs.idx_) :
                difference_type(rhs.idx_ - idx_);
        }

    protected:
        Container* data_;
        size_t idx_;
    };


    template<typename, typename, typename, typename, typename>
    class const_stable_iterator;


    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer,
             typename Distance  = typename Container::difference_type>
    class stable_iterator :
        public stable_iterator_base<stable_iterator<Container, ValueType, Reference, Pointer, Distance>,
                                    Container, ValueType, Reference, Pointer, Distance>
    {
    public:
        using _my_base = stable_iterator_base<stable_iterator, Container, ValueType, Reference, Pointer, Distance>;
        using _my_base::_my_base;

        friend class const_stable_iterator<Container, ValueType, Reference, Pointer, Distance>;
    };


    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer,
             typename Distance  = typename Container::difference_type>
    class const_stable_iterator :
        public stable_iterator_base<const_stable_iterator<Container, ValueType, Reference, Pointer, Distance>,
                                    const Container, const ValueType, const Reference, const Pointer, Distance>
    {
    public:
        using _my_base = stable_iterator_base<const_stable_iterator, const Container, const ValueType, const Reference, const Pointer, Distance>;
        using _my_base::_my_base;

        /* implicit */ const_stable_iterator(stable_iterator<Container, ValueType, Reference, Pointer, Distance> it) noexcept :
            _my_base(*it.data_, it.idx_)
        {}
    };


    template<typename Container>
    inline auto stable_begin(Container& container) noexcept
    {
        return stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    inline auto stable_end(Container& container) noexcept
    {
        return stable_iterator<Container>(container, container.size());
    }

    template<typename Container>
    inline auto stable_begin(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    inline auto stable_end(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, container.size());
    }

    template<typename Container>
    inline auto stable_cbegin(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    inline auto stable_cend(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, container.size());
    }


    /* iota iterator */

    template<std::integral T = std::size_t, typename Difference = std::ptrdiff_t>
    class iota_iterator : public random_access_iterator_interface<iota_iterator<T>, Difference>
    {
    public:
        using _my_base = random_access_iterator_interface<iota_iterator, Difference>;

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using value_type = T;
        using reference  = T;
        using pointer    = T;

        iota_iterator() noexcept :
            value_(T{})
        {}

        iota_iterator(const T& val) noexcept :
            value_(val)
        {}

        pointer operator->() const = delete;

        reference operator*() const
        {
            return value_;
        }

        bool operator==(const iota_iterator& rhs) const
        {
            return value_ == rhs.value_;
        }

        bool operator<(const iota_iterator& rhs) const
        {
            return value_ < rhs.value_;
        }

        iota_iterator& increment() noexcept
        {
            GA_ASSERT(value_ != std::numeric_limits<T>::max(), "Can't increment iterator with max value.");

            ++value_;
            return *this;
        }

        iota_iterator& decrement() noexcept
        {
            GA_ASSERT(value_ != std::numeric_limits<T>::min(), "Can't decrement iterator with min value.");

            --value_;
            return *this;
        }

        iota_iterator& operator+=(difference_type n)
        {
            // TODO ar conversions for max/min?
            GA_ASSERT(n > 0 ? (std::numeric_limits<T>::max() - n) >= value_ : true, "Can't increment iterator past its max value.");
            GA_ASSERT(n < 0 ? (std::numeric_limits<T>::min() - n) <= value_ : true, "Can't decrement iterator past its min value.");

            value_ += n;
            return *this;
        }

        difference_type operator-(const iota_iterator& rhs) const
        {
            // TODO bounds checks
            return difference_type(value_) - difference_type(rhs.value_);
        }

    private:
        T value_;
    };

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP