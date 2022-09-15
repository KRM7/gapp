/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ITERATORS_HPP
#define GA_UTILITY_ITERATORS_HPP

#include <iterator>

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
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(this); }
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
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(this); }
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

        bool operator!=(const Derived& rhs) const
        {
            return !(*this == rhs);
        }

    private:
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(this); }
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
        constexpr Derived& derived() noexcept { return *static_cast<Derived*>(this); }
        constexpr const Derived& derived() const noexcept { return *static_cast<const Derived*>(this); }
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(this); }
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
        constexpr Derived& derived() noexcept { return static_cast<Derived&>(this); }
        constexpr const Derived& derived() const noexcept { return static_cast<const Derived&>(this); }
    };

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP