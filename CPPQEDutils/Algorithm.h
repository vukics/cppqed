// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Generic algorithms not found in either STL or Boost}
#ifndef   CPPQEDCORE_UTILS_ALGORITHM_H_INCLUDED
#define   CPPQEDCORE_UTILS_ALGORITHM_H_INCLUDED

#include <boost/range/numeric.hpp>

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



namespace cppqedutils {


/// Fills a container by output iterator with concatenated values taken subsequently from the input sequences. \note Concatenation can be expressed as accumulation
template<typename SeqOfSeqs, typename Out_Iterator>
const Out_Iterator
concatenateViaIterator(const SeqOfSeqs& sOs, ///<[in] the sequence containing the input sequences
                       Out_Iterator out)
{
  return boost::accumulate(sOs,out, [](auto iter, const auto& t){return std::copy(t.begin(),t.end(),iter);});
}


/// Fills a container of the necessary size with concatenated values taken subsequently from the input sequences
template<typename SeqOfSeqs, typename Out>
const Out&
concatenate(const SeqOfSeqs& sOs, ///<[in] the sequence containing the input sequences
            Out& out)
{
  concatenateViaIterator(sOs,out.begin());
  return out;
}


/// Fills an *empty* (default-constructed) container with concatenated values taken subsequently from the input sequences
template<typename Out, typename SeqOfSeqs>
const Out
concatenateGrow(const SeqOfSeqs& sOs)
{
  Out empty;
  concatenateViaIterator(sOs,std::back_inserter(empty));
  return empty;
}


} // cppqedutils


#endif // CPPQEDCORE_UTILS_ALGORITHM_H_INCLUDED
