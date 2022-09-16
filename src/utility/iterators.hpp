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

        auto operator[](difference_type n) const { return *(*this + n); }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    template<typename Derived,
             typename ValueType,
             typename Reference,
             typename Pointer>
    class stable_iterator_base : public random_access_iterator_interface<Derived>
    {
    public:

        using _my_base = random_access_iterator_interface<Derived>;

        using iterator_category = std::random_access_iterator_tag;

        using value_type = ValueType;
        using reference  = Reference;
        using pointer    = Pointer;
        using typename _my_base::difference_type;

        reference operator*() const
        {
            assert(derived().data_ && derived().data_->size() > derived().idx_);
            return derived().data_[derived().idx_];
        }

        friend bool operator==(Derived lhs, Derived rhs)
        {
            assert(lhs.data_ == rhs.data_);
            return (lhs.data_ == rhs.data_) && (lhs.idx_ == rhs.idx_);
        }

        friend bool operator<(Derived lhs, Derived rhs)
        {
            assert(lhs.data_ == rhs.data_);
            return lhs.idx_ < rhs.idx_;
        }

        Derived& operator++()
        {
            assert(derived().data_);
            derived().idx_++;
            return derived();
        }

        Derived& operator--()
        {
            assert(derived().data_ && derived().idx_ != 0);
            derived().idx_--;
            return derived();
        }

        Derived& operator+=(difference_type n)
        {
            assert(derived().data_); assert(n < 0 ? (derived().idx_ >= -n) : true);
            derived().idx_ += n;
            return derived();
        }

        friend difference_type operator-(Derived lhs, Derived rhs)
        {
            assert(lhs.data_ && lhs.data_ == rhs.data_);
            return difference_type(lhs.idx_) - difference_type(rhs.idx_);
        }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(*this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(*this); }
    };


    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer>
    class stable_iterator : public stable_iterator_base<stable_iterator<Container, ValueType, Reference, Pointer>,
                                                        ValueType, Reference, Pointer>
    {
    public:
        using _my_base = stable_iterator_base<stable_iterator, ValueType, Reference, Pointer>;

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using typename _my_base::value_type;
        using typename _my_base::reference;
        using typename _my_base::pointer;

        stable_iterator() noexcept :
            _my_base(), data_(nullptr), idx_(0)
        {}

        stable_iterator(Container& container, size_t idx) noexcept :
            _my_base(), data_(&container), idx_(idx)
        {}

    private:
        Container* data_;
        size_t idx_;

        friend class _my_base;
    };


    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer>
    class const_stable_iterator : public stable_iterator_base<const_stable_iterator<Container, ValueType, Reference, Pointer>,
                                                              const ValueType, const Reference, const Pointer>
    {
    public:
        using _my_base = stable_iterator_base<const_stable_iterator, const ValueType, const Reference, const Pointer>;

        using typename _my_base::iterator_category;
        using typename _my_base::difference_type;
        using typename _my_base::value_type;
        using typename _my_base::reference;
        using typename _my_base::pointer;

        const_stable_iterator() noexcept :
            data_(nullptr), idx_(0)
        {}

        const_stable_iterator(const Container& container, size_t idx) noexcept :
            data_(&container), idx_(idx)
        {}

        const_stable_iterator(stable_iterator<Container, ValueType, Reference, Pointer> it) noexcept :
            data_(it.data_), idx_(it.idx)
        {}

    private:
        const Container* data_;
        size_t idx_;

        friend class _my_base;
    };


    /*template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer>
    class stable_iterator : public random_access_iterator_interface<stable_iterator<Container>>
    {
    public:
        using iterator_category = std::random_access_iterator_tag;

        using typename random_access_iterator_interface<stable_iterator>::difference_type;
        using value_type = ValueType;
        using reference  = Reference;
        using pointer    = Pointer;

        stable_iterator() noexcept :
            data_(nullptr), idx_(0)
        {}

        stable_iterator(Container& container, size_t idx) noexcept :
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

    
    template<typename Container,
             typename ValueType = typename Container::value_type,
             typename Reference = typename Container::reference,
             typename Pointer   = typename Container::pointer>
    class const_stable_iterator : public random_access_iterator_interface<const_stable_iterator<Container>>
    {
    public:
        using iterator_category = std::random_access_iterator_tag;

        using typename random_access_iterator_interface<const_stable_iterator>::difference_type;
        using value_type = const ValueType;
        using reference  = const Reference;
        using pointer    = const Pointer;

        const_stable_iterator() noexcept :
            data_(nullptr), idx_(0)
        {}

        const_stable_iterator(const Container& container, size_t idx) noexcept :
            data_(&container), idx_(idx)
        {}

        const_stable_iterator(stable_iterator<Container, ValueType, Reference, Pointer> it) noexcept :
            data_(it.data_), idx_(it.idx)
        {}

        reference operator*() const
        {
            assert(data_ && data_->size() > idx_);
            return data_[idx_];
        }

        friend bool operator==(const_stable_iterator lhs, const_stable_iterator rhs)
        {
            assert(lhs.data_ == rhs.data_);
            return (lhs.data_ == rhs.data_) && (lhs.idx_ == rhs.idx_);
        }

        friend bool operator<(const_stable_iterator lhs, const_stable_iterator rhs)
        {
            assert(lhs.data_ == rhs.data_);
            return lhs.idx_ < rhs.idx_;
        }

        const_stable_iterator& operator++()
        {
            assert(data_);
            ++idx_;
            return *this;
        }

        const_stable_iterator& operator--()
        {
            assert(data_ && idx_ != 0);
            --idx_;
            return *this;
        }

        const_stable_iterator& operator+=(difference_type n)
        {
            assert(data_); assert(n < 0 ? (idx_ >= -n) : true);
            idx_ += n;
            return *this;
        }

        friend difference_type operator-(const_stable_iterator lhs, const_stable_iterator rhs)
        {
            assert(lhs.data_ && lhs.data_ == rhs.data_);
            return difference_type(lhs.idx_) - difference_type(rhs.idx_);
        }
        
    private:
        const Container* data_;
        size_t idx_;
    };*/


    template<typename Container>
    auto stable_begin(Container& container) noexcept
    {
        return stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    auto stable_end(Container& container) noexcept
    {
        return stable_iterator<Container>(container, container.size());
    }

    template<typename Container>
    auto stable_begin(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, 0);
    }

    template<typename Container>
    auto stable_end(const Container& container) noexcept
    {
        return const_stable_iterator<Container>(container, container.size());
    }

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP