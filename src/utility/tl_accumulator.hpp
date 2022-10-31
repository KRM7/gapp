/* Copyright (c) 2022 Krisztián Rugási. Subject to the MIT License. */

#ifndef GA_UTILITY_TL_ACCUMULATOR_HPP
#define GA_UTILITY_TL_ACCUMULATOR_HPP

#include <vector>
#include <algorithm>
#include <functional>
#include <mutex>
#include <cstddef>

namespace genetic_algorithm::detail
{
    template<typename T, size_t Instance = 0, typename ReduceOp = std::plus<T>>
    class tl_vector_accumulator
    {
    public:
        static decltype(auto) at(size_t i) noexcept
        {
            return tl_instance_.data_[i];
        }

        /* Note: Calling this while concurrently writing to a reference returned by at() is undefined behaviour. */
        static std::vector<T> collect()
        {
            std::scoped_lock _(lock_);

            std::vector<T> sum_vec(size_);

            for (tl_vector_accumulator* tl_vec : tl_vector_list_)
            {
                std::transform(tl_vec->data_.begin(), tl_vec->data_.end(), sum_vec.begin(), sum_vec.begin(), ReduceOp{});
            }
            std::transform(accumulator_vector_.begin(), accumulator_vector_.end(), sum_vec.begin(), sum_vec.begin(), ReduceOp{});

            return sum_vec;
        }

        /* Note: Calling this while concurrently writing to a reference returned by at() is undefined behaviour. */
        static void reset(size_t size, const T& initial_value = T{})
        {
            std::scoped_lock _(lock_);

            size_ = size;
            for (tl_vector_accumulator* vec : tl_vector_list_)
            {
                vec->data_.resize(size_);
                std::fill(vec->data_.begin(), vec->data_.end(), initial_value);
            }
            accumulator_vector_.resize(size_);
            std::fill(accumulator_vector_.begin(), accumulator_vector_.end(), T{});
        }

        static size_t size() noexcept
        {
            std::scoped_lock _(lock_);
            return size_;
        }

    private:
        tl_vector_accumulator() /* on thread init */
        {
            std::scoped_lock _(lock_);
            tl_vector_list_.push_back(this);
        }

        ~tl_vector_accumulator() /* on thread exit */
        {
            std::scoped_lock _(lock_);
            std::transform(data_.begin(), data_.end(), accumulator_vector_.begin(), accumulator_vector_.begin(), ReduceOp{});
            std::erase(tl_vector_list_, this);
        }

        inline static size_t size_;
        inline static std::vector<tl_vector_accumulator*> tl_vector_list_;
        inline static std::vector<T> accumulator_vector_;
        inline static std::mutex lock_;

        static thread_local tl_vector_accumulator tl_instance_;
        std::vector<T> data_ = std::vector<T>(size_);
    };

    template<typename T, size_t I, typename R>
    inline thread_local tl_vector_accumulator<T, I, R> tl_vector_accumulator<T, I, R>::tl_instance_;

} // namespace genetic_algorithm::detail

#endif // !GA_UTILITY_TL_ACCUMULATOR_HPP