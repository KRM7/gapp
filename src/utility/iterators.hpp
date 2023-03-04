/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ITERATORS_HPP
#define GA_UTILITY_ITERATORS_HPP

#include "utility.hpp"
#include "type_traits.hpp"
#include <iterator>
#include <type_traits>
#include <concepts>
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
        constexpr auto cbegin() const noexcept { return derived().begin(); };
        constexpr auto cend() const noexcept { return derived().end(); };

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
        constexpr auto rbegin()        { return std::make_reverse_iterator(derived().end()); }
        constexpr auto rbegin() const  { return std::make_reverse_iterator(derived().end()); }
        constexpr auto crbegin() const { return std::make_reverse_iterator(derived().cend()); }
        constexpr auto rend()          { return std::make_reverse_iterator(derived().begin()); }
        constexpr auto rend() const    { return std::make_reverse_iterator(derived().begin()); }
        constexpr auto crend() const   { return std::make_reverse_iterator(derived().cbegin()); }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
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

        constexpr Derived& operator++()
        {
            return derived().increment();
        }

        constexpr Derived operator++(int)
        {
            Derived old_value(derived());
            derived().increment();
            return old_value;
        }

        constexpr auto operator->() const
        {
            static_assert(std::is_lvalue_reference_v<dereference_t<Derived>>);
            return &*derived();
        }

        friend constexpr bool operator!=(const Derived& lhs, const Derived& rhs)
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

        constexpr Derived& operator--()
        {
            return derived().decrement();
        }

        constexpr Derived operator--(int)
        {
            Derived old_value(derived());
            derived().decrement();
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

        friend constexpr bool operator>(const Derived& lhs, const Derived& rhs)  { return rhs < lhs; }
        friend constexpr bool operator<=(const Derived& lhs, const Derived& rhs) { return !(rhs < lhs); }
        friend constexpr bool operator>=(const Derived& lhs, const Derived& rhs) { return !(lhs < rhs); }

        friend constexpr Derived operator+(Derived iter, difference_type n) { return iter += n; }
        friend constexpr Derived operator+(difference_type n, Derived iter) { return iter += n; }

        friend constexpr Derived& operator-=(Derived& lhs, difference_type n) { return lhs += -n; }
        friend constexpr Derived operator-(Derived iter, difference_type n) { return iter -= n; }

        constexpr decltype(auto) operator[](difference_type n) const { return *(derived() + n); }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
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
    class contiguous_iterator_interface : public random_access_iterator_interface<Derived, Distance>
    {
    public:
        using iterator_category = std::contiguous_iterator_tag;
    };


    /* Iterators for random-access containers that aren't invalidated on reallocations. */

    template<typename Derived,
             typename Container,
             typename ValueType,
             typename Reference,
             typename Pointer,
             typename Distance = std::ptrdiff_t>
    class stable_iterator_base : public random_access_iterator_interface<Derived, Distance>
    {
    public:
        using my_base_ = random_access_iterator_interface<Derived>;

        using typename my_base_::iterator_category;
        using typename my_base_::difference_type;
        using value_type = ValueType;
        using reference  = Reference;
        using pointer    = Pointer;

        constexpr stable_iterator_base() noexcept :
            data_(nullptr), idx_(0)
        {}

        constexpr stable_iterator_base(Container& container, size_t idx) noexcept :
            data_(&container), idx_(idx)
        {
            GA_ASSERT(data_->size() >= idx_, "Iterator can't refer to element past the end of the range.");
        }

        constexpr reference operator*() const
        {
            GA_ASSERT(data_ != nullptr, "Can't dereference value initialized iterator.");
            GA_ASSERT(data_->size() > idx_, "Can't dereference past-the-end iterator.");

            return (*data_)[idx_];
        }

        friend constexpr bool operator==(const Derived& lhs, const Derived& rhs)
        {
            GA_ASSERT(lhs.data_ == rhs.data_, "Can't compare iterators of different ranges.");
            GA_ASSERT(lhs.data_ == nullptr || lhs.data_->size() >= lhs.idx_, "Can't compare invalid iterator.");
            GA_ASSERT(rhs.data_ == nullptr || rhs.data_->size() >= rhs.idx_, "Can't compare invalid iterator.");

            return lhs.idx_ == rhs.idx_;  /* Value-initialized iterators will have the same idx. */
        }

        friend constexpr bool operator<(const Derived& lhs, const Derived& rhs)
        {
            GA_ASSERT(lhs.data_ == rhs.data_, "Can't compare iterators of different ranges.");
            GA_ASSERT(lhs.data_ == nullptr || lhs.data_->size() >= lhs.idx_, "Can't compare invalid iterator.");
            GA_ASSERT(rhs.data_ == nullptr || rhs.data_->size() >= rhs.idx_, "Can't compare invalid iterator.");

            return lhs.idx_ < rhs.idx_;   /* Value-initialized iterators will have the same idx. */
        }

        constexpr Derived& increment()
        {
            GA_ASSERT(data_ != nullptr, "Can't increment value initialized iterator.");
            GA_ASSERT(idx_ != data_->size(), "Can't increment past-the-end iterator.");

            ++idx_;
            return static_cast<Derived&>(*this);
        }

        constexpr Derived& decrement()
        {
            GA_ASSERT(data_ != nullptr, "Can't decrement value initialized iterator.");
            GA_ASSERT(idx_ != 0, "Can't decremenet the begin iterator.");

            --idx_;
            return static_cast<Derived&>(*this);
        }

        constexpr Derived& operator+=(difference_type n)
        {
            GA_ASSERT(data_ != nullptr, "Can't offset value initialized iterator.");
            GA_ASSERT(n < 0 ? idx_ >= size_t(-n) : true, "Can't move iterator to before the start of the range.");
            GA_ASSERT(n > 0 ? idx_ <= (data_->size() - n) : true, "Can't move iterator past the end of the range.");

            idx_ += n;
            return static_cast<Derived&>(*this);
        }

        friend constexpr difference_type operator-(const Derived& lhs, const Derived& rhs)
        {
            GA_ASSERT(lhs.data_ && rhs.data_, "Can't get the distance of value initialized iterators.");
            GA_ASSERT(lhs.data_ == rhs.data_, "Can't get the distance of iterators of different ranges.");
            GA_ASSERT(lhs.data_->size() >= lhs.idx_, "Invalid lhs iterator.");
            GA_ASSERT(rhs.data_->size() >= rhs.idx_, "Invalid rhs iterator.");

            const size_t distance = lhs.idx_ >= rhs.idx_ ?
                lhs.idx_ - rhs.idx_ :
                rhs.idx_ - lhs.idx_;

            GA_ASSERT(distance <= size_t(std::numeric_limits<difference_type>::max()), "Can't represent the result of the operation as difference_type.");

            return difference_type(distance);
        }

    protected:
        Container* data_;
        size_t idx_;
    };


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
        using my_base_ = stable_iterator_base<stable_iterator, Container, ValueType, Reference, Pointer, Distance>;
        using my_base_::my_base_;

        template<typename, typename, typename, typename, typename, typename>
        friend class const_stable_iterator;
    };


    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::const_reference,
             typename Pointer   = typename Container::const_pointer,
             typename Distance  = typename Container::difference_type,
             typename Iterator  = stable_iterator<Container, ValueType, typename Container::reference, typename Container::pointer, Distance>>
    class const_stable_iterator :
        public stable_iterator_base<const_stable_iterator<Container, ValueType, Reference, Pointer, Distance>,
                                    const Container, ValueType, Reference, Pointer, Distance>
    {
    public:
        using my_base_ = stable_iterator_base<const_stable_iterator, const Container, ValueType, Reference, Pointer, Distance>;
        using my_base_::my_base_;

        /* implicit */ constexpr const_stable_iterator(Iterator it) noexcept
        {
            this->data_ = it.data_;
            this->idx_  = it.idx_;
        }
    };


    template<typename Container>
    constexpr auto stable_begin(Container& container) noexcept
    {
        return stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    constexpr auto stable_end(Container& container) noexcept
    {
        return stable_iterator<Container>(container, container.size());
    }

    template<typename Container>
    constexpr auto stable_begin(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    constexpr auto stable_end(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, container.size());
    }

    template<typename Container>
    constexpr auto stable_cbegin(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    constexpr auto stable_cend(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, container.size());
    }


    /* iota iterator */

    template<std::integral T = std::size_t>
    class iota_iterator : public random_access_iterator_interface<iota_iterator<T>, std::make_signed_t<T>>
    {
    public:
        using my_base_ = random_access_iterator_interface<iota_iterator, std::make_signed_t<T>>;

        using typename my_base_::iterator_category;
        using typename my_base_::difference_type;
        using value_type = T;
        using reference  = T;
        using pointer    = T;

        constexpr iota_iterator() noexcept :
            value_()
        {}

        explicit constexpr iota_iterator(T val) noexcept :
            value_(val)
        {}

        pointer operator->() const = delete;

        constexpr reference operator*() const noexcept
        {
            return value_;
        }

        friend constexpr bool operator==(const iota_iterator& lhs, const iota_iterator& rhs) noexcept
        {
            return lhs.value_ == rhs.value_;
        }

        friend constexpr bool operator<(const iota_iterator& lhs, const iota_iterator& rhs) noexcept
        {
            return lhs.value_ < rhs.value_;
        }

        constexpr iota_iterator& increment()
        {
            GA_ASSERT(value_ != std::numeric_limits<T>::max(), "Can't increment iterator with max value.");

            ++value_;
            return *this;
        }

        constexpr iota_iterator& decrement()
        {
            GA_ASSERT(value_ != std::numeric_limits<T>::min(), "Can't decrement iterator with min value.");

            --value_;
            return *this;
        }

        constexpr iota_iterator& operator+=(difference_type n)
        {
            GA_ASSERT(n > 0 ? (std::numeric_limits<T>::max() - n) >= value_ : true, "Can't increment iterator past its max value.");
            GA_ASSERT(n < 0 ? size_t(-n) <= (value_ - std::numeric_limits<T>::min()) : true, "Can't decrement iterator past its min value.");

            value_ += n;
            return *this;
        }

        friend constexpr difference_type operator-(iota_iterator lhs, iota_iterator rhs)
        {
            const T distance = lhs.value_ >= rhs.value_ ?
                lhs.value_ - rhs.value_ :
                rhs.value_ - lhs.value_;

            GA_ASSERT(distance <= T(std::numeric_limits<difference_type>::max()), "Can't represent the result of the operation as difference_type.");

            return difference_type(distance);
        }

    private:
        T value_;
    };

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP