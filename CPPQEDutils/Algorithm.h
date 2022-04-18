// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Range algorithms not yet in STL in C++20}
#pragma once

#include <iterator>
#include <ranges>


namespace cppqedutils::ranges {

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
  

} // cppqedutils::ranges

