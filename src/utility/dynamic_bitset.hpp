/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_DYNAMIC_BITSET_HPP
#define GAPP_UTILITY_DYNAMIC_BITSET_HPP

#include "small_vector.hpp"
#include "bit.hpp"
#include "utility.hpp"
#include <memory>
#include <bit>
#include <limits>
#include <cstdint>
#include <cstddef>

// NOLINTBEGIN(*bool-conversion, *assignment, *assignment-signature, *operator, *ref-data-members)

namespace gapp::detail
{
    class dynamic_bitset
    {
    public:
        using value_type      = bool;
        using block_type      = std::size_t;
        using allocator_type  = std::allocator<block_type>;
        using const_reference = bool;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        static constexpr size_type block_size = std::numeric_limits<block_type>::digits;

        class reference
        {
        public:
            constexpr /* implicit */ operator bool() const noexcept
            {
                return static_cast<bool>(block_ & mask_);
            }

            constexpr const reference& operator=(bool value) const noexcept
            {
                value ? set() : clear();
                return *this;
            }

            constexpr const reference& operator=(const reference& rhs) const noexcept
            {
                return *this = static_cast<bool>(rhs);
            }

            constexpr void set() const noexcept   { block_ |= mask_; }
            constexpr void clear() const noexcept { block_ &= ~mask_; }
            constexpr void flip() const noexcept  { block_ ^= mask_; }
            
            constexpr void operator&() = delete;

        private:
            constexpr reference(dynamic_bitset& bitset, size_type idx) noexcept :
                block_(bitset.blocks_[idx / block_size]),
                mask_(block_type(1) << (idx % block_size))
            {
                GAPP_ASSERT(idx < bitset.size());
            }

            friend dynamic_bitset;

            block_type& block_;
            block_type mask_;
        };

        constexpr dynamic_bitset() = default;

        constexpr dynamic_bitset(size_type size) :
            blocks_(size / block_size + bool(size % block_size)),
            size_(size)
        {}

        constexpr dynamic_bitset(size_type size, bool value) :
            blocks_(size / block_size + bool(size % block_size), block_of<block_type>(value)),
            size_(size)
        {}

        dynamic_bitset(const dynamic_bitset&)            = default;
        dynamic_bitset(dynamic_bitset&&)                 = default;
        dynamic_bitset& operator=(const dynamic_bitset&) = default;
        dynamic_bitset& operator=(dynamic_bitset&&)      = default;

        constexpr reference operator[](size_type idx) noexcept
        {
            GAPP_ASSERT(idx < size_);

            return { *this, idx };
        }

        constexpr const_reference operator[](size_type idx) const noexcept
        {
            GAPP_ASSERT(idx < size_);

            const size_type block  = idx / block_size;
            const size_type offset = idx % block_size;

            return static_cast<bool>(blocks_[block] & (block_type(1) << offset));
        }

        constexpr bool empty() const noexcept
        {
            return size_ == 0;
        }

        constexpr size_type size() const noexcept
        {
            return size_;
        }

        constexpr void clear() noexcept
        {
            blocks_.clear();
            size_ = 0;
        }

        constexpr void resize(size_type new_size, bool value = false)
        {
            blocks_.resize(new_size / block_size + bool(new_size % block_size), block_of<block_type>(value));

            if (new_size > size_)
            {
                blocks_.back() &= partial_block_mask();
                blocks_.back() |= ~partial_block_mask() & block_of<block_type>(value);
            }

            size_ = new_size;
        }

        constexpr void fill(bool value) noexcept
        {
            for (block_type& block : blocks_)
            {
                block = block_of<block_type>(value);
            }
        }

        constexpr size_type find_first(bool value) const noexcept
        {
            return value ? find_first_one() : find_first_zero();
        }

        constexpr size_type popcount() const noexcept
        {
            size_type count = 0;

            for (size_type i = 0; i < full_block_count(); i++)
            {
                count += std::popcount(blocks_[i]);
            }

            return count + std::popcount(partial_block());
        }

        constexpr bool any_set() const noexcept
        {
            for (size_type i = 0; i < full_block_count(); i++)
            {
                if (blocks_[i] != zeros<block_type>) return true;
            }

            return partial_block() != zeros<block_type>;
        }

        constexpr bool all_set() const noexcept
        {
            for (size_type i = 0; i < full_block_count(); i++)
            {
                if (blocks_[i] != ones<block_type>) return false;
            }

            return size_type(std::popcount(partial_block())) == partial_block_size();
        }

        constexpr bool none_set() const noexcept
        {
            return !any_set();
        }

        constexpr friend dynamic_bitset operator~(const dynamic_bitset& bitset) noexcept
        {
            dynamic_bitset complement(bitset.size());

            for (size_type i = 0; i < bitset.blocks_.size(); i++)
            {
                complement.blocks_[i] = ~bitset.blocks_[i];
            }

            return complement;
        }

        constexpr friend bool operator==(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            if (lhs.size() != rhs.size()) return false;

            for (std::size_t i = 0; i < lhs.full_block_count(); i++)
            {
                if (lhs.blocks_[i] != rhs.blocks_[i]) return false;
            }

            return lhs.partial_block() == rhs.partial_block();
        }

    private:
        small_vector<block_type, 4> blocks_;
        size_type size_ = 0;

        constexpr size_type find_first_one() const noexcept
        {
            for (size_type i = 0; i < blocks_.size(); i++)
            {
                if (blocks_[i] == zeros<block_type>) continue;
                return std::min(size_, i * block_size + std::countr_zero(blocks_[i]));
            }
            return size_;
        }

        constexpr size_type find_first_zero() const noexcept
        {
            for (size_type i = 0; i < blocks_.size(); i++)
            {
                if (blocks_[i] == ones<block_type>) continue;
                return std::min(size_, i * block_size + std::countr_one(blocks_[i]));
            }
            return size_;
        }

        constexpr size_type full_block_count() const noexcept
        {
            return size_ / block_size;
        }

        constexpr size_type partial_block_size() const noexcept
        {
            return size_ % block_size;
        }

        constexpr block_type partial_block_mask() const noexcept
        {
            return mask_right_n<block_type>(partial_block_size());
        }

        constexpr block_type partial_block() const noexcept
        {
            return blocks_.back() & partial_block_mask();
        }
    };

} // namespace gapp::detail

// NOLINTEND(*bool-conversion, *assignment, *assignment-signature, *operator, *ref-data-members)

#endif // !GAPP_UTILITY_DYNAMIC_BITSET_HPP
