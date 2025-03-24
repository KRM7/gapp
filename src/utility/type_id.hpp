/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_TYPE_ID_HPP
#define GAPP_UTILITY_TYPE_ID_HPP

#include <bit>
#include <cstddef>

namespace gapp::detail
{
    template<typename T>
    inline size_t type_id() noexcept
    {
        constexpr static unsigned char dummy = 0;
        return std::bit_cast<size_t>(&dummy);
    }

} // namespace gapp::detail

#endif // !GAPP_UTILITY_TYPE_ID_HPP
