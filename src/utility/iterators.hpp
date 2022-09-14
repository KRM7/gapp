/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_ITERATORS_HPP
#define GA_UTILITY_ITERATORS_HPP

#include <iterator>

namespace genetic_algorithm::detail
{
    template<typename Derived>
    class reverse_iterator_interface
    {
    public:
        /* Requires begin(), end(), cbegin(), cend() to be implemented in Derived. */
        auto rbegin() noexcept        { return std::make_reverse_iterator(derived().end()); }
        auto rbegin() const noexcept  { return std::make_reverse_iterator(derived().end()); }
        auto crbegin() const noexcept { return std::make_reverse_iterator(derived().cend()); }
        auto rend() noexcept          { return std::make_reverse_iterator(derived().begin()); }
        auto rend() const noexcept    { return std::make_reverse_iterator(derived().begin()); }
        auto crend() const noexcept   { return std::make_reverse_iterator(derived().cbegin()); }

    private:
        constexpr Derived& derived() noexcept { return *static_cast<Derived*>(this); }
        constexpr const Derived& derived() const noexcept { return *static_cast<const Derived*>(this); }
    };
    

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_ITERATORS_HPP