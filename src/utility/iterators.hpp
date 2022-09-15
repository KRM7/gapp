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

    private:
        constexpr Derived& derived() noexcept { return *static_cast<Derived*>(this); }
        constexpr const Derived& derived() const noexcept { return *static_cast<const Derived*>(this); }
    };
    

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP