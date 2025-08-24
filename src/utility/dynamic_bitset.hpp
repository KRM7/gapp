/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_DYNAMIC_BITSET_HPP
#define GAPP_UTILITY_DYNAMIC_BITSET_HPP

#include "small_vector.hpp"
#include "iterators.hpp"
#include "bit.hpp"
#include "hash.hpp"
#include "utility.hpp"
#include <algorithm>
#include <iterator>
#include <memory>
#include <bit>
#include <compare>
#include <limits>
#include <cstdint>
#include <cstddef>

// NOLINTBEGIN(*bool-conversion, *assignment, *assignment-signature, *operator, *ref-data-members)

namespace gapp
{
    class dynamic_bitset : public detail::iterator_interface<dynamic_bitset>
    {
    public:
        using value_type      = bool;
        using block_type      = std::size_t;
        using allocator_type  = std::allocator<block_type>;
        using const_reference = bool;
        using pointer         = void;
        using const_pointer   = void;
        using size_type       = std::size_t;
        using difference_type = std::ptrdiff_t;

        class reference;

        using iterator       = detail::stable_iterator<dynamic_bitset>;
        using const_iterator = detail::const_stable_iterator<dynamic_bitset>;

        using reverse_iterator       = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

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

            constexpr const reference& operator&=(bool value) const noexcept
            {
                if (!value) clear();
                return *this;
            }

            constexpr const reference& operator|=(bool value) const noexcept
            {
                if (value) set();
                return *this;
            }

            constexpr const reference& operator^=(bool value) const noexcept
            {
                if (value) flip();
                return *this;
            }

            constexpr const reference& operator-=(bool value) const noexcept
            {
                if (value) clear();
                return *this;
            }

            constexpr void set() const noexcept   { block_ |= mask_; }
            constexpr void clear() const noexcept { block_ &= ~mask_; }
            constexpr void flip() const noexcept  { block_ ^= mask_; }

            constexpr void swap(const reference& other) const noexcept
            {
                if (*this != other)
                {
                    this->flip();
                    other.flip();
                }
            }

            constexpr friend void swap(const reference& lhs, const reference& rhs) noexcept
            {
                lhs.swap(rhs);
            }
            
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

        constexpr explicit dynamic_bitset(size_type size) :
            blocks_(size / block_size + bool(size % block_size)),
            size_(size)
        {}

        constexpr dynamic_bitset(size_type size, bool value) :
            blocks_(size / block_size + bool(size % block_size), detail::block_of<block_type>(value)),
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

        constexpr iterator begin() noexcept { return { this, 0 }; }
        constexpr iterator end() noexcept { return { this, size_ }; }

        constexpr const_iterator begin() const noexcept { return { this, 0 }; }
        constexpr const_iterator end() const noexcept { return { this, size_ }; }

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
            blocks_.resize((new_size + block_size) / block_size, detail::block_of<block_type>(value));

            if (new_size > size_)
            {
                blocks_.back() &= partial_block_mask();
                blocks_.back() |= ~partial_block_mask() & detail::block_of<block_type>(value);
            }

            size_ = new_size;
        }

        constexpr void push_back(bool value)
        {
            if (size_ % block_size == 0) blocks_.push_back(block_type{ 0 });

            (*this)[size_++] = value;
        }

        constexpr void fill(bool value) noexcept
        {
            for (block_type& block : blocks_)
            {
                block = detail::block_of<block_type>(value);
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
                if (blocks_[i] != detail::zeros<block_type>) return true;
            }

            return partial_block() != detail::zeros<block_type>;
        }

        constexpr bool all_set() const noexcept
        {
            for (size_type i = 0; i < full_block_count(); i++)
            {
                if (blocks_[i] != detail::ones<block_type>) return false;
            }

            return size_type(std::popcount(partial_block())) == partial_block_size();
        }

        constexpr bool none_set() const noexcept
        {
            return !any_set();
        }

        constexpr std::span<block_type> blocks() noexcept
        {
            return blocks_;
        }

        constexpr std::span<const block_type> blocks() const noexcept
        {
            return blocks_;
        }

        constexpr void swap(dynamic_bitset& other) noexcept
        {
            blocks_.swap(other.blocks_);
            std::swap(size_, other.size_);
        }

        constexpr friend void swap(dynamic_bitset& lhs, dynamic_bitset& rhs) noexcept
        {
            lhs.swap(rhs);
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

        constexpr friend auto operator<=>(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            const size_type block_count = std::min(lhs.full_block_count(), rhs.full_block_count());

            for (size_t i = 0; i < block_count; i++)
            {
                const block_type lhs_block = lhs.blocks_[i];
                const block_type rhs_block = rhs.blocks_[i];

                if (lhs_block != rhs_block) return lhs_block <=> rhs_block;
            }

            if (lhs.size() < rhs.size())
            {
                const block_type lhs_block = lhs.partial_block();
                const block_type rhs_block = rhs.blocks_[block_count] & lhs.partial_block_mask();

                if (lhs_block != rhs_block) return lhs_block <=> rhs_block;

                return std::strong_ordering::less;
            }

            if (lhs.size() > rhs.size())
            {
                const block_type lhs_block = lhs.blocks_[block_count] & rhs.partial_block_mask();
                const block_type rhs_block = rhs.partial_block();

                if (lhs_block != rhs_block) return lhs_block <=> rhs_block;

                return std::strong_ordering::greater;
            }

            GAPP_ASSERT(lhs.size() == rhs.size());

            return lhs.partial_block() <=> rhs.partial_block();
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

        constexpr void flip() noexcept
        {
            for (block_type& block : blocks_) { block = ~block; }
        }

        constexpr friend dynamic_bitset operator&(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            dynamic_bitset result(lhs.size());

            for (size_type i = 0; i < result.blocks_.size(); i++)
            {
                result.blocks_[i] = lhs.blocks_[i] & rhs.blocks_[i];
            }

            return result;
        }

        constexpr friend dynamic_bitset& operator&=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            for (size_type i = 0; i < lhs.blocks_.size(); i++)
            {
                lhs.blocks_[i] &= rhs.blocks_[i];
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator|(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            dynamic_bitset result(lhs.size());

            for (size_type i = 0; i < result.blocks_.size(); i++)
            {
                result.blocks_[i] = lhs.blocks_[i] | rhs.blocks_[i];
            }

            return result;
        }

        constexpr friend dynamic_bitset& operator|=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            for (size_type i = 0; i < lhs.blocks_.size(); i++)
            {
                lhs.blocks_[i] |= rhs.blocks_[i];
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator^(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            dynamic_bitset result(lhs.size());

            for (size_type i = 0; i < result.blocks_.size(); i++)
            {
                result.blocks_[i] = lhs.blocks_[i] ^ rhs.blocks_[i];
            }

            return result;
        }

        constexpr friend dynamic_bitset& operator^=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            for (size_type i = 0; i < lhs.blocks_.size(); i++)
            {
                lhs.blocks_[i] ^= rhs.blocks_[i];
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator-(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            dynamic_bitset result(lhs.size());

            for (size_type i = 0; i < result.blocks_.size(); i++)
            {
                result.blocks_[i] = lhs.blocks_[i] & ~rhs.blocks_[i];
            }

            return result;
        }

        constexpr friend dynamic_bitset& operator-=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            for (size_type i = 0; i < lhs.blocks_.size(); i++)
            {
                lhs.blocks_[i] &= ~rhs.blocks_[i];
            }

            return lhs;
        }

    private:
        small_vector<block_type, 4> blocks_;
        size_type size_ = 0;

        constexpr size_type find_first_one() const noexcept
        {
            for (size_type i = 0; i < blocks_.size(); i++)
            {
                if (blocks_[i] == detail::zeros<block_type>) continue;
                return std::min(size_, i * block_size + std::countr_zero(blocks_[i]));
            }
            return size_;
        }

        constexpr size_type find_first_zero() const noexcept
        {
            for (size_type i = 0; i < blocks_.size(); i++)
            {
                if (blocks_[i] == detail::ones<block_type>) continue;
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
            return detail::mask_right_n<block_type>(partial_block_size());
        }

        constexpr block_type partial_block() const noexcept
        {
            return !empty() ? blocks_.back() & partial_block_mask() : block_type{ 0 };
        }
    };

} // namespace gapp

namespace std
{
    template<>
    struct hash<gapp::dynamic_bitset>
    {
        std::size_t operator()(const gapp::dynamic_bitset& bitset) const noexcept
        {
            return gapp::detail::hash_range(bitset.blocks());
        }
    };

} // namespace std

// NOLINTEND(*bool-conversion, *assignment, *assignment-signature, *operator, *ref-data-members)

#endif // !GAPP_UTILITY_DYNAMIC_BITSET_HPP
