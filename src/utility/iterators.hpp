/* Copyright (c) 2022 Kriszti�n Rug�si. Subject to the MIT License. */

#ifndef GA_UTILITY_ITERATORS_HPP
#define GA_UTILITY_ITERATORS_HPP

#include "utility.hpp"
#include <iterator>
#include <type_traits>
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

        auto operator->() const requires (std::is_lvalue_reference_v<decltype(*this->derived())>)
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
            GA_ASSERT(data_, "Can't dereference value initialized iterator.");
            GA_ASSERT(data_->size() > idx_, "Can't dereference past-the-end iterator.");

            return (*data_)[idx_];
        }

        bool operator==(Derived rhs) const
        {
            GA_ASSERT(data_ == rhs.data_, "Can't compare iterators of different ranges.");

            return (idx_ == rhs.idx_) || (data_ == nullptr && rhs.data_ == nullptr);
        }

        bool operator<(Derived rhs) const
        {
            GA_ASSERT(data_ == rhs.data_, "Can't compare iterators of different ranges.");

            return idx_ < rhs.idx_;
        }

        Derived& increment() noexcept
        {
            GA_ASSERT(data_, "Can't increment value initialized iterator.");

            idx_++;
            return derived();
        }

        Derived& decrement() noexcept
        {
            GA_ASSERT(data_, "Can't decrement value initialized iterator.");

            idx_--;
            return derived();
        }

        Derived& operator+=(difference_type n)
        {
            GA_ASSERT(derived().data_, "Can't offset value initialized iterator.");
            GA_ASSERT(n < 0 ? (difference_type(idx_) >= -n) : true, "Can't move iterator to before the start of the range.");
            GA_ASSERT(n > 0 ? (idx_ + n) <= data_->size() : true, "Can't move iterator past the end of the range.");

            idx_ += n;
            return derived();
        }

        difference_type operator-(Derived rhs) const
        {
            GA_ASSERT(data_ && rhs.data_, "Can't get the difference of value initialized iterators.");
            GA_ASSERT(data_ == rhs.data_, "Can't get the distance of iterators of different ranges");

            return difference_type(idx_) - difference_type(rhs.idx_);
        }

    protected:

        Container* data_;
        size_t idx_;

    private:

        Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
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

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using typename _my_base::value_type;
        using typename _my_base::reference;
        using typename _my_base::pointer;

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

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using typename _my_base::value_type;
        using typename _my_base::reference;
        using typename _my_base::pointer;

        using _my_base::_my_base;

        /* implicit */ const_stable_iterator(stable_iterator<Container, ValueType, Reference, Pointer, Distance> it) noexcept :
            _my_base()
        {
            this->data_ = it.data_;
            this->idx_ = it.idx_;
        }
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

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP