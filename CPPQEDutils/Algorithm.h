// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Range algorithms not yet in STL in C++20}
#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>


namespace std_ext::ranges {

using namespace std;

/// ranges::fold and co. based on [this proposal](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2020/p2214r0.html#stdaccumulate-rangesfold)
template <class F, class T, class U>
concept foldable =
  regular_invocable<F&, T, U> &&
  convertible_to<invoke_result_t<F&, T, U>, T>;

template <class F, class T, class I>
concept indirectly_binary_foldable =
  indirectly_readable<I> &&
  copy_constructible<F> &&
  foldable<F, T, iter_value_t<I>> &&
  foldable<F, T, iter_reference_t<I>> &&
  foldable<F, T, iter_common_reference_t<I>>;


template <::std::ranges::input_range R, movable T, class Proj = identity, indirectly_binary_foldable<T, projected<::std::ranges::iterator_t<R>, Proj>> BinaryOperation>
constexpr T fold(R&& r, T init, BinaryOperation op, Proj proj = {}) {
  ::std::ranges::iterator_t<R> b = begin(r);
  ::std::ranges::sentinel_t<R> e = end(r);
  for (; b != e; ++b) {
      init = op(std::move(init), proj(*b));
  }
  return init;
}
  

} // std_ext::ranges


namespace cppqedutils {


/// std::array`s of different size cannot be std::views::join`ed because the size is part of the type
/** this solution based on C++17 fold expressions comes from [here](http://stackoverflow.com/a/42774523/1171157) */
template <typename Type, std::size_t... sizes>
constexpr auto concatenate(const std::array<Type, sizes>&... arrays)
{
  std::array<Type, (sizes + ...)> result;
  std::size_t index{};

  ((std::copy_n(arrays.begin(), sizes, result.begin() + index), index += sizes), ...);

  return result;
}


} // cppqedutils
