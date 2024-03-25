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
                return static_cast<bool>(*block_ & mask_);
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

            constexpr void set() const noexcept { *block_ |= mask_; }
            constexpr void clear() const noexcept { *block_ &= ~mask_; }
            constexpr void flip() const noexcept { *block_ ^= mask_; }
            
            constexpr void operator&() = delete;

        private:
            constexpr reference(dynamic_bitset& bitset, size_type idx) noexcept :
                block_(std::addressof(bitset.blocks_[idx / block_size])),
                mask_(block_type{ 1 } << (idx % block_size))
            {}

            friend dynamic_bitset;

            block_type* block_;
            block_type mask_;
        };

        constexpr dynamic_bitset() = default;

        constexpr dynamic_bitset(size_type size) :
            blocks_(size / block_size + bool(size % block_size)),
            size_(size)
        {}

        constexpr dynamic_bitset(size_type size, bool value) :
            blocks_(size / block_size + bool(size % block_size), to_mask<block_type>(value)),
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

            return static_cast<bool>(blocks_[block] & (block_type{ 1 } << offset));
        }

        constexpr bool empty() const noexcept { return size_ == 0; }
        constexpr size_type size() const noexcept { return size_; }

        constexpr void clear() noexcept
        {
            blocks_.clear();
            size_ = 0;
        }

        constexpr void resize(size_type size, bool value = false)
        {
            blocks_.resize(size / block_size + bool(size % block_size), to_mask<block_type>(value));

            if (size > size_)
            {
                const size_type count = block_size - size_ % block_size;
                const block_type mask = ones_left_n<block_type>(count);

                blocks_.back() &= ~mask;
                blocks_.back() |= mask & to_mask<block_type>(value);
            }

            size_ = size;
        }

        constexpr void fill(bool value) noexcept
        {
            for (block_type& block : blocks_)
            {
                block = to_mask<block_type>(value);
            }
        }

        constexpr size_type find_first(bool value) const noexcept
        {
            return value ? find_first_one() : find_first_zero();
        }

        constexpr size_type popcount() const noexcept
        {
            size_type count = 0;

            for (size_type i = 0; i < size_ / block_size; i++)
            {
                count += std::popcount(blocks_[i]);
            }

            if (size_type rem_size = size_ % block_size; rem_size)
            {
                count += std::popcount(blocks_.back() & ones_right_n<block_type>(rem_size));
            }

            return count;
        }

        constexpr bool any_set() const noexcept
        {
            for (size_type i = 0; i < size_ / block_size; i++)
            {
                if (blocks_[i] != zeros<block_type>) return true;
            }

            if (size_type rem_size = size_ % block_size; rem_size)
            {
                const block_type rem_block = blocks_.back() & ones_right_n<block_type>(rem_size);
                if (rem_block != zeros<block_type>) return true;
            }

            return false;
        }

        constexpr bool all_set() const noexcept
        {
            for (size_type i = 0; i < size_ / block_size; i++)
            {
                if (blocks_[i] != ones<block_type>) return false;
            }

            if (size_type rem_size = size_ % block_size; rem_size)
            {
                const block_type rem_block = blocks_.back() | ones_left_n<block_type>(block_size - rem_size);
                if (rem_block != ones<block_type>) return false;
            }

            return true;
        }

        constexpr bool none_set() const noexcept
        {
            return !any_set();
        }

        constexpr friend dynamic_bitset operator~(const dynamic_bitset& bitset) noexcept
        {
            dynamic_bitset result(bitset.size());

            for (size_type i = 0; i < bitset.blocks_.size(); i++)
            {
                result.blocks_[i] = ~bitset.blocks_[i];
            }

            return result;
        }

        constexpr friend bool operator==(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            if (lhs.size_ != rhs.size_) return false;

            for (std::size_t i = 0; i < lhs.size_ / lhs.block_size; i++)
            {
                if (lhs.blocks_[i] != rhs.blocks_[i]) return false;
            }

            if (std::size_t rem_size = lhs.size_ % lhs.block_size; rem_size)
            {
                auto rem_mask = ones_right_n<typename dynamic_bitset::block_type>(rem_size);
                auto lhs_rem_block = lhs.blocks_.back() & rem_mask;
                auto rhs_rem_block = rhs.blocks_.back() & rem_mask;

                if (lhs_rem_block != rhs_rem_block) return false;
            }

            return true;
        }

    private:
        small_vector<block_type> blocks_;
        size_type size_ = 0;

        constexpr size_type find_first_one() const noexcept
        {
            for (size_type i = 0; i < blocks_.size(); i++)
            {
                if (blocks_[i] != zeros<block_type>) return std::min(size_, i * block_size + std::countr_zero(blocks_[i]));
            }
            return size_;
        }

        constexpr size_type find_first_zero() const noexcept
        {
            for (size_type i = 0; i < blocks_.size(); i++)
            {
                if (blocks_[i] != ones<block_type>) return std::min(size_, i * block_size + std::countr_one(blocks_[i]));
            }
            return size_;
        }
    };

} // namespace gapp::detail

#endif // !GAPP_UTILITY_DYNAMIC_BITSET_HPP
