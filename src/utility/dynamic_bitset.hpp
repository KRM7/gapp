/* Copyright (c) 2024 Krisztián Rugási. Subject to the MIT License. */

#ifndef GAPP_UTILITY_DYNAMIC_BITSET_HPP
#define GAPP_UTILITY_DYNAMIC_BITSET_HPP

#include "small_vector.hpp"
#include "iterators.hpp"
#include "functional.hpp"
#include "bit.hpp"
#include "hash.hpp"
#include "utility.hpp"
#include <span>
#include <algorithm>
#include <iterator>
#include <memory>
#include <bit>
#include <compare>
#include <limits>
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

        constexpr explicit dynamic_bitset(size_type size, bool value = false) :
            blocks_((size + block_size - 1) / block_size, detail::block_of<block_type>(value)),
            size_(size)
        {
            if (partial_block_size()) blocks_.back() &= partial_block_mask();

            GAPP_ASSERT(unused_bits_are_zero());
        }

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
            if (value && new_size > size_ && partial_block_size())
            {
                blocks_.back() |= ~partial_block_mask();
            }

            blocks_.resize((new_size + block_size - 1) / block_size, detail::block_of<block_type>(value));
            size_ = new_size;

            if (partial_block_size())
            {
                blocks_.back() &= partial_block_mask();
            }

            GAPP_ASSERT(unused_bits_are_zero());
        }

        constexpr void push_back(bool value)
        {
            if (!partial_block_size()) blocks_.push_back(0);

            (*this)[size_++] = value;
        }

        constexpr void fill(bool value) noexcept
        {
            auto first = blocks_.begin();
            auto last  = blocks_.begin() + full_block_count();

            for ( ; first != last; first++)
            {
                *first = detail::block_of<block_type>(value);
            }

            if (partial_block_size())
            {
                blocks_.back() = detail::block_of<block_type>(value) & partial_block_mask();
            }

            GAPP_ASSERT(unused_bits_are_zero());
        }

        constexpr size_type find_first(bool value) const noexcept
        {
            return value ? find_first_one() : find_first_zero();
        }

        constexpr size_type popcount() const noexcept
        {
            GAPP_ASSERT(unused_bits_are_zero());

            size_type count = 0;

            for (block_type block : blocks_)
            {
                count += std::popcount(block);
            }

            return count;
        }

        constexpr bool any_set() const noexcept
        {
            GAPP_ASSERT(unused_bits_are_zero());

            return std::any_of(blocks_.begin(), blocks_.end(), detail::not_equal_to(block_type{ 0 }));
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
            GAPP_ASSERT(unused_bits_are_zero());

            return blocks_;
        }

        constexpr std::span<const block_type> blocks() const noexcept
        {
            GAPP_ASSERT(unused_bits_are_zero());

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
            GAPP_ASSERT(lhs.unused_bits_are_zero());
            GAPP_ASSERT(rhs.unused_bits_are_zero());

            if (lhs.size() != rhs.size()) return false;

            const size_type block_count = lhs.blocks_.size();

            auto lhs_block = lhs.blocks_.begin();
            auto rhs_block = rhs.blocks_.begin();

            for (size_type i = 0; i < block_count; i++)
            {
                if (*lhs_block++ != *rhs_block++) return false;
            }

            return true;
        }

        constexpr friend auto operator<=>(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            const size_type block_count = std::min(lhs.full_block_count(), rhs.full_block_count());

            for (size_type i = 0; i < block_count; i++)
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
            dynamic_bitset result(bitset);
            return result.flip();
        }

        constexpr dynamic_bitset& flip() noexcept
        {
            auto first = blocks_.begin();
            auto last  = blocks_.begin() + full_block_count();

            for ( ; first != last; first++)
            {
                *first = ~*first;
            }

            if (partial_block_size())
            {
                blocks_.back() = ~blocks_.back() & partial_block_mask();
            }

            return *this;
        }

        constexpr friend dynamic_bitset operator&(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            dynamic_bitset result(lhs);
            return result &= rhs;
        }

        constexpr friend dynamic_bitset& operator&=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            GAPP_ASSERT(lhs.unused_bits_are_zero());
            GAPP_ASSERT(rhs.unused_bits_are_zero());

            const size_type block_count = lhs.blocks_.size();

            auto lhs_block = lhs.blocks_.begin();
            auto rhs_block = rhs.blocks_.begin();

            for (size_type i = 0; i < block_count; i++)
            {
                *lhs_block++ &= *rhs_block++;
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator|(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            dynamic_bitset result(lhs);
            return result |= rhs;
        }

        constexpr friend dynamic_bitset& operator|=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            GAPP_ASSERT(lhs.unused_bits_are_zero());
            GAPP_ASSERT(rhs.unused_bits_are_zero());

            const size_type block_count = lhs.blocks_.size();

            auto lhs_block = lhs.blocks_.begin();
            auto rhs_block = rhs.blocks_.begin();

            for (size_type i = 0; i < block_count; i++)
            {
                *lhs_block++ |= *rhs_block++;
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator^(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            dynamic_bitset result(lhs);
            return result ^= rhs;
        }

        constexpr friend dynamic_bitset& operator^=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            GAPP_ASSERT(lhs.unused_bits_are_zero());
            GAPP_ASSERT(rhs.unused_bits_are_zero());

            const size_type block_count = lhs.blocks_.size();

            auto lhs_block = lhs.blocks_.begin();
            auto rhs_block = rhs.blocks_.begin();

            for (size_type i = 0; i < block_count; i++)
            {
                *lhs_block++ ^= *rhs_block++;
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator-(const dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            dynamic_bitset result(lhs);
            return result -= rhs;
        }

        constexpr friend dynamic_bitset& operator-=(dynamic_bitset& lhs, const dynamic_bitset& rhs) noexcept
        {
            GAPP_ASSERT(lhs.size() == rhs.size());

            GAPP_ASSERT(lhs.unused_bits_are_zero());
            GAPP_ASSERT(rhs.unused_bits_are_zero());

            const size_type block_count = lhs.blocks_.size();

            auto lhs_block = lhs.blocks_.begin();
            auto rhs_block = rhs.blocks_.begin();

            for (size_type i = 0; i < block_count; i++)
            {
                *lhs_block++ &= ~*rhs_block++;
            }

            return lhs;
        }

        constexpr friend dynamic_bitset operator>>(const dynamic_bitset& lhs, size_type n) noexcept
        {
            dynamic_bitset result(lhs);
            return result >>= n;
        }

        constexpr friend dynamic_bitset& operator>>=(dynamic_bitset& lhs, size_type n) noexcept
        {
            GAPP_ASSERT(n <= lhs.size());
            GAPP_ASSERT(lhs.unused_bits_are_zero());

            if (n == 0) return lhs;

            const size_type block_disp = n / block_size;
            const size_type shift_disp = n % block_size;

            auto& blocks = lhs.blocks_;

            if (!shift_disp)
            {
                std::move(blocks.begin() + block_disp, blocks.end(), blocks.begin());
                std::fill(blocks.end() - block_disp, blocks.end(), block_type{ 0 });
                return lhs;
            }

            for (size_type i = 0; i < blocks.size() - block_disp - 1; i++)
            {
                const size_type lo_idx = i + block_disp;
                const size_type hi_idx = i + block_disp + 1;
                blocks[i] = (blocks[lo_idx] >> shift_disp) | (blocks[hi_idx] << (block_size - shift_disp));
            }

            blocks[blocks.size() - block_disp - 1] = (blocks.back() >> shift_disp);

            std::fill(blocks.end() - block_disp, blocks.end(), block_type{ 0 });
            return lhs;
        }

        constexpr friend dynamic_bitset operator<<(const dynamic_bitset& lhs, size_type n) noexcept
        {
            dynamic_bitset result(lhs);
            return result <<= n;
        }

        constexpr friend dynamic_bitset& operator<<=(dynamic_bitset& lhs, size_t n) noexcept
        {
            GAPP_ASSERT(n <= lhs.size());
            GAPP_ASSERT(lhs.unused_bits_are_zero());

            if (n == 0) return lhs;

            const size_type block_disp = n / block_size;
            const size_type shift_disp = n % block_size;

            auto& blocks = lhs.blocks_;

            if (!shift_disp)
            {
                std::move_backward(blocks.begin(), blocks.end() - block_disp, blocks.begin() + block_disp);
                std::fill(blocks.begin(), blocks.begin() + block_disp, block_type{ 0 });
                return lhs;
            }

            for (size_type i = blocks.size() - block_disp - 1; i > 0; i--)
            {
                const size_type lo_idx = i - 1;
                const size_type hi_idx = i;
                blocks[i + block_disp] = (blocks[hi_idx] << shift_disp) | (blocks[lo_idx] >> (block_size - shift_disp));
            }

            blocks[block_disp] = (blocks[0] << shift_disp);
            blocks.back() &= lhs.partial_block_mask();

            std::fill(blocks.begin(), blocks.begin() + block_disp, block_type{ 0 });
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
            GAPP_ASSERT(unused_bits_are_zero());

            return partial_block_size() ? blocks_.back() : block_type{0};
        }

        constexpr bool unused_bits_are_zero() const noexcept
        {
            return !partial_block_size() || !(blocks_.back() & ~partial_block_mask());
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
