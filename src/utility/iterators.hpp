/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_ITERATORS_HPP
#define GAPP_UTILITY_ITERATORS_HPP

#include "utility.hpp"
#include "type_traits.hpp"
#include <iterator>
#include <ranges>
#include <type_traits>
#include <concepts>
#include <compare>
#include <limits>
#include <memory>
#include <stdexcept>
#include <cstddef>

namespace gapp::detail
{
    /*
    * The following should be implemented in Derived:
    *   begin(), end() with const overloads
    */
    template<typename Derived>
    class iterator_interface
    {
    public:
        constexpr auto cbegin() const noexcept //requires(std::input_iterator<detail::iterator_t<Derived>>)
        {
            return derived().begin();
        }

        constexpr auto cend() const noexcept //requires(std::input_iterator<detail::iterator_t<Derived>>)
        {
            return derived().end();
        }

        constexpr auto rbegin() //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return std::make_reverse_iterator(derived().end());
        }

        constexpr auto rbegin() const //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return std::make_reverse_iterator(derived().end());
        }

        constexpr auto crbegin() const //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return std::make_reverse_iterator(derived().end());
        }

        constexpr auto rend() //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return std::make_reverse_iterator(derived().begin());
        }

        constexpr auto rend() const //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return std::make_reverse_iterator(derived().begin());
        }

        constexpr auto crend() const //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return std::make_reverse_iterator(derived().begin());
        }

        friend auto operator<=>(const iterator_interface&, const iterator_interface&) = default;

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /*
    * The following should be implemented in Derived:
    *   begin(), end() with const overloads
    */
    template<typename Derived>
    class container_interface : public iterator_interface<Derived>
    {
    public:
        constexpr size_t size() const noexcept //requires(std::forward_iterator<detail::iterator_t<Derived>>)
        {
            return std::distance(derived().begin(), derived().end());
        }

        constexpr bool empty() const noexcept //requires(std::input_iterator<detail::iterator_t<Derived>>)
        {
            return derived().begin() == derived().end();
        }

        constexpr decltype(auto) front() //requires(std::input_iterator<detail::iterator_t<Derived>>)
        {
            return *derived().begin();
        }

        constexpr decltype(auto) front() const //requires(std::input_iterator<detail::iterator_t<Derived>>)
        {
            return *derived().begin();
        }

        constexpr decltype(auto) back() //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return *std::prev(derived().end());
        }

        constexpr decltype(auto) back() const //requires(std::bidirectional_iterator<detail::iterator_t<Derived>>)
        {
            return *std::prev(derived().end());
        }

        constexpr decltype(auto) operator[](size_t idx) //requires(std::random_access_iterator<detail::iterator_t<Derived>>)
        {
            GAPP_ASSERT(idx < size(), "Index out of bounds.");

            return *std::next(derived().begin(), idx);
        }

        constexpr decltype(auto) operator[](size_t idx) const //requires(std::random_access_iterator<detail::iterator_t<Derived>>)
        {
            GAPP_ASSERT(idx < size(), "Index out of bounds.");

            return *std::next(derived().begin(), idx);
        }

        constexpr decltype(auto) at(size_t idx) //requires(std::random_access_iterator<detail::iterator_t<Derived>>)
        {
            if (idx >= size()) throw std::out_of_range{ "Index out of bounds." };

            return *std::next(derived().begin(), idx);
        }

        constexpr decltype(auto) at(size_t idx) const //requires(std::random_access_iterator<detail::iterator_t<Derived>>)
        {
            if (idx >= size()) throw std::out_of_range{ "Index out of bounds." };

            return *std::next(derived().begin(), idx);
        }

        constexpr auto data() noexcept //requires(std::contiguous_iterator<detail::iterator_t<Derived>>)
        {
            return std::to_address(derived().begin());
        }

        constexpr auto data() const noexcept //requires(std::contiguous_iterator<detail::iterator_t<Derived>>)
        {
            return std::to_address(derived().begin());
        }

        friend auto operator<=>(const container_interface&, const container_interface&) = default;

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    /* Helper class for the iterators' operator-> for when it doesn't return an lvalue reference. */
    template<typename T>
    class proxy_ptr
    {
    public:
        constexpr proxy_ptr(T data) noexcept(std::is_nothrow_move_constructible_v<T>) :
            data_(std::move(data))
        {}

        constexpr T* operator->() noexcept { return std::addressof(data_); }

    private:
        T data_;
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
            std::ignore = derived().increment();
            return old_value;
        }

        constexpr auto operator->() const
        {
            if constexpr (std::is_lvalue_reference_v<dereference_t<Derived>>) return std::addressof(*derived());
            else return detail::proxy_ptr{ *derived() };
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
            std::ignore = derived().decrement();
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
    *   operator+=(n)
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
    *   operator+=(n)
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

        constexpr stable_iterator_base(Container* container, size_t idx) noexcept :
            data_(container), idx_(idx)
        {}

        constexpr stable_iterator_base(Container& container, size_t idx) noexcept :
            data_(std::addressof(container)), idx_(idx)
        {
            GAPP_ASSERT(data_->size() >= idx_, "Iterator can't refer to element past the end of the range.");
        }

        constexpr reference operator*() const
        {
            GAPP_ASSERT(data_ != nullptr, "Can't dereference value initialized iterator.");
            GAPP_ASSERT(data_->size() > idx_, "Can't dereference past-the-end iterator.");

            return (*data_)[idx_];
        }

        friend constexpr bool operator==(const Derived& lhs, const Derived& rhs)
        {
            GAPP_ASSERT(lhs.data_ == rhs.data_, "Can't compare iterators of different ranges.");
            GAPP_ASSERT(lhs.data_ == nullptr || lhs.data_->size() >= lhs.idx_, "Can't compare invalid lhs iterator.");
            GAPP_ASSERT(rhs.data_ == nullptr || rhs.data_->size() >= rhs.idx_, "Can't compare invalid rhs iterator.");

            return lhs.idx_ == rhs.idx_;  /* Value-initialized iterators will have the same idx. */
        }

        friend constexpr bool operator<(const Derived& lhs, const Derived& rhs)
        {
            GAPP_ASSERT(lhs.data_ == rhs.data_, "Can't compare iterators of different ranges.");
            GAPP_ASSERT(lhs.data_ == nullptr || lhs.data_->size() >= lhs.idx_, "Can't compare invalid lhs iterator.");
            GAPP_ASSERT(rhs.data_ == nullptr || rhs.data_->size() >= rhs.idx_, "Can't compare invalid rhs iterator.");

            return lhs.idx_ < rhs.idx_;   /* Value-initialized iterators will have the same idx. */
        }

        constexpr Derived& increment()
        {
            GAPP_ASSERT(data_ != nullptr, "Can't increment value initialized iterator.");
            GAPP_ASSERT(idx_ != data_->size(), "Can't increment past-the-end iterator.");

            ++idx_;
            return static_cast<Derived&>(*this);
        }

        constexpr Derived& decrement()
        {
            GAPP_ASSERT(data_ != nullptr, "Can't decrement value initialized iterator.");
            GAPP_ASSERT(idx_ != 0, "Can't decrement the begin iterator.");

            --idx_;
            return static_cast<Derived&>(*this);
        }

        constexpr Derived& operator+=(difference_type n)
        {
            GAPP_ASSERT(data_ != nullptr, "Can't offset value initialized iterator.");
            GAPP_ASSERT(n < 0 ? idx_ >= size_t(-n) : true, "Can't move iterator to before the start of the range.");
            GAPP_ASSERT(n > 0 ? idx_ <= (data_->size() - n) : true, "Can't move iterator past the end of the range.");

            idx_ += n;
            return static_cast<Derived&>(*this);
        }

        friend constexpr difference_type operator-(const Derived& lhs, const Derived& rhs)
        {
            GAPP_ASSERT(lhs.data_ && rhs.data_, "Can't get the distance of value initialized iterators.");
            GAPP_ASSERT(lhs.data_ == rhs.data_, "Can't get the distance of iterators of different ranges.");
            GAPP_ASSERT(lhs.data_->size() >= lhs.idx_, "Invalid lhs iterator.");
            GAPP_ASSERT(rhs.data_->size() >= rhs.idx_, "Invalid rhs iterator.");

            const size_t distance = lhs.idx_ >= rhs.idx_ ?
                lhs.idx_ - rhs.idx_ :
                rhs.idx_ - lhs.idx_;

            GAPP_ASSERT(distance <= size_t(std::numeric_limits<difference_type>::max()), "Can't represent the result of the operation as difference_type.");

            return difference_type(distance);
        }

    protected:
        Container* data_;
        size_t idx_;
    };


    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::const_reference,
             typename Pointer   = typename Container::const_pointer,
             typename Distance  = typename Container::difference_type>
    class const_stable_iterator :
        public stable_iterator_base<const_stable_iterator<Container, ValueType, Reference, Pointer, Distance>,
                                    const Container, ValueType, Reference, Pointer, Distance>
    {
    public:
        using my_base_ = stable_iterator_base<const_stable_iterator, const Container, ValueType, Reference, Pointer, Distance>;
        using my_base_::my_base_;
    };

    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer,
             typename Distance  = typename Container::difference_type,
             typename CIterator = const_stable_iterator<Container, ValueType, typename Container::const_reference, typename Container::const_pointer, Distance>>
    class stable_iterator :
        public stable_iterator_base<stable_iterator<Container, ValueType, Reference, Pointer, Distance>,
                                    Container, ValueType, Reference, Pointer, Distance>
    {
    public:
        using my_base_ = stable_iterator_base<stable_iterator, Container, ValueType, Reference, Pointer, Distance>;
        using my_base_::my_base_;
        using const_iterator = CIterator;

        constexpr /* implicit */ operator const_iterator() const noexcept { return { this->data_, this->idx_ }; }
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
        using pointer    = detail::proxy_ptr<T>;

        constexpr iota_iterator() noexcept :
            value_()
        {}

        explicit constexpr iota_iterator(T val) noexcept :
            value_(val)
        {}

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
            GAPP_ASSERT(value_ != std::numeric_limits<T>::max(), "Can't increment iterator with max value.");

            ++value_;
            return *this;
        }

        constexpr iota_iterator& decrement()
        {
            GAPP_ASSERT(value_ != std::numeric_limits<T>::min(), "Can't decrement iterator with min value.");

            --value_;
            return *this;
        }

        constexpr iota_iterator& operator+=(difference_type n)
        {
            GAPP_ASSERT(n > 0 ? (std::numeric_limits<T>::max() - n) >= value_ : true, "Can't increment iterator past its max value.");
            GAPP_ASSERT(n < 0 ? size_t(-n) <= (value_ - std::numeric_limits<T>::min()) : true, "Can't decrement iterator past its min value.");

            value_ += n;
            return *this;
        }

        friend constexpr difference_type operator-(iota_iterator lhs, iota_iterator rhs)
        {
            return (lhs.value_ >= rhs.value_) ? difference_type(lhs.value_ - rhs.value_) : -difference_type(rhs.value_ - lhs.value_);
        }

    private:
        T value_;
    };


    /* base iterator */

    template<typename Base>
    class base_iterator : public random_access_iterator_interface<base_iterator<Base>>
    {
    public:
        using my_base_ = random_access_iterator_interface<base_iterator<Base>>;

        using typename my_base_::iterator_category;
        using typename my_base_::difference_type;

        using value_type = Base;
        using reference  = Base&;
        using pointer    = Base*;

        constexpr base_iterator() noexcept = default;

        template<std::derived_from<Base> Derived>
        explicit constexpr base_iterator(Derived* data) noexcept :
            ptr_(reinterpret_cast<internal_pointer>(static_cast<Base*>(data))),
            step_(sizeof(Derived))
        {}

        constexpr reference operator*() const noexcept
        {
            return *std::launder(reinterpret_cast<pointer>(ptr_));
        }

        friend constexpr bool operator==(const base_iterator& lhs, const base_iterator& rhs) noexcept
        {
            return lhs.ptr_ == rhs.ptr_;
        }

        friend constexpr bool operator<(const base_iterator& lhs, const base_iterator& rhs) noexcept
        {
            return lhs.ptr_ < rhs.ptr_;
        }

        constexpr base_iterator& increment() noexcept
        {
            ptr_ += step_;
            return *this;
        }

        constexpr base_iterator& decrement() noexcept
        {
            ptr_ -= step_;
            return *this;
        }

        constexpr base_iterator& operator+=(difference_type n) noexcept
        {
            ptr_ += n * step_;
            return *this;
        }

        friend constexpr difference_type operator-(base_iterator lhs, base_iterator rhs) noexcept
        {
            GAPP_ASSERT(lhs.step_ == rhs.step_);
            return (lhs.ptr_ - rhs.ptr_) / lhs.step_;
        }

    private:
        using internal_pointer = detail::copy_cv_t<Base, unsigned char>*;

        internal_pointer ptr_ = nullptr;
        std::size_t step_ = 0;
    };


    template<typename Base, std::ranges::contiguous_range Container>
    constexpr auto base_begin(Container& container) noexcept
    {
        if (container.empty()) return base_iterator<Base>();

        return base_iterator<Base>(std::to_address(container.begin()));
    }

    template<typename Base, std::ranges::contiguous_range Container>
    constexpr auto base_end(Container& container) noexcept
    {
        return base_begin<Base>(container) + container.size();
    }

    template<typename Base>
    struct base_view : public container_interface<base_view<Base>>
    {
        template<std::ranges::contiguous_range Range>
        requires(std::is_lvalue_reference_v<Range>)
        constexpr /* implicit */ base_view(Range&& r) noexcept :
            first_(base_begin<Base>(r)),
            last_(base_end<Base>(r))
        {}

        constexpr auto begin() const noexcept { return first_; }
        constexpr auto end() const noexcept { return last_; }

    private:
        base_iterator<Base> first_;
        base_iterator<Base> last_;
    };


    template<std::input_iterator Iter>
    constexpr auto make_move_iterator_if_noexcept(Iter it) noexcept
    {
        if constexpr (std::is_nothrow_move_constructible_v<std::iter_value_t<Iter>>)
            return std::make_move_iterator(it);
        else return it;
    }

} // namespace gapp::detail

#endif // !GAPP_UTILITY_ITERATORS_HPP
