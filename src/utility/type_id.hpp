/* Copyright (c) 2023 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_TYPE_ID_HPP
#define GA_UTILITY_TYPE_ID_HPP

#include <cstddef>

namespace genetic_algorithm::detail
{
    class TypeIdGen final
    {
        static size_t next() noexcept
        {
            static size_t id = 0;
            return id++;
        }

        template<typename T>
        friend struct TypeId;
    };

    template<typename T>
    struct TypeId final
    {
        inline static const size_t value = TypeIdGen::next();
    };

    template<typename T>
    inline const size_t type_id = TypeId<T>::value;

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_TYPE_ID_HPP
