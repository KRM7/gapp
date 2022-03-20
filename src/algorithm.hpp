#ifndef GA_ALGORITHM_HPP
#define GA_ALGORITHM_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <type_traits>
#include <concepts>
#include <cstddef>
#include <cassert>

namespace genetic_algorithm::detail
{
    template<typename T>
    inline constexpr auto lforward(T&& t) noexcept
    {
        return std::conditional_t<std::is_lvalue_reference_v<T>,
                                  std::reference_wrapper<std::remove_reference_t<T>>,
                                  T>
               { std::forward<T>(t) };
    }


    template<typename F>
    inline constexpr auto compose(F&& f) noexcept
    {
        return [f = lforward(f)] <typename... Args>
        (Args&&... args) requires std::invocable<F, Args...>
        {
            return std::invoke(f, std::forward<Args>(args)...);
        };
    }

    template<typename F, typename... Fs>
    inline constexpr auto compose(F&& f, Fs&&... fs) noexcept
    {
        return [f = lforward(f), ...fs = lforward(fs)] <typename... Args>
        (Args&&... args) requires std::invocable<F, Args...>
        {
            return compose(fs...)(std::invoke(f, std::forward<Args>(args)...));
        };
    }

    template<typename Container, typename F>
    requires std::invocable<F, const typename Container::value_type>
    inline auto map(const Container& cont, F&& f)
    {
        using result_t = std::vector<std::invoke_result_t<F, const typename Container::value_type>>;

        result_t result;
        result.reserve(cont.size());

        std::transform(std::begin(cont), std::end(cont), std::back_inserter(result),
        [f = lforward(f)](const Container::value_type& elem)
        {
            return std::invoke(f, elem);
        });

        return result;
    }

    template<std::random_access_iterator Iter,
             typename Comp = std::less<typename std::iterator_traits<Iter>::value_type>>
    requires std::invocable<Comp, const typename std::iterator_traits<Iter>::value_type,
                                  const typename std::iterator_traits<Iter>::value_type>
    std::vector<size_t> argsort(Iter first, Iter last, Comp&& comp = std::less<typename std::iterator_traits<Iter>::value_type>{})
    {
        assert(std::distance(first, last) >= 0);

        std::vector<size_t> indices(std::distance(first, last));
        std::iota(indices.begin(), indices.end(), size_t{ 0 });

        std::sort(indices.begin(), indices.begin(),
        [first, last, comp = lforward(comp)](size_t lidx, size_t ridx)
        {
            return std::invoke(comp, *std::next(first, lidx), *std::next(first, ridx));
        });

        return indices;
    }

    // argmin(first, last, comp) -> [size_t idx, val_t val];
}

#endif // !GA_ALGORITHM_HPP