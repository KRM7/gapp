/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_TYPE_ID_HPP
#define GA_UTILITY_TYPE_ID_HPP

#include <bit>
#include <cstddef>

namespace gapp::detail
{
    template<typename T>
    struct TypeIdHelper { constexpr static int var = 0; };

    template<typename T>
    inline size_t type_id() noexcept
    {
        return std::bit_cast<size_t>(&TypeIdHelper<T>::var);
    }

} // namespace gapp::detail

#endif // !GA_UTILITY_TYPE_ID_HPP
