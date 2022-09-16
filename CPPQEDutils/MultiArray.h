// Copyright András Vukics 2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Algorithm.h"
#include "TMP_Tools.h"

#include <boost/range/combine.hpp>

#include <boost/json.hpp>

#include <array>
#include <concepts>
#include <functional>
#include <span>
#include <stdexcept>
#include <vector>


// TODO: RetainedAxes can be a non-type template parameter in C++20!!! And maybe retainedAxes can be a std::index_sequence, which can be converted to hana::tuple when needed?

namespace cppqedutils {


/// should be passed by value
template <size_t RANK>
using Extents=std::array<size_t,RANK>;


template <size_t RANK>
auto& incrementMultiIndex(Extents<RANK>& idx, Extents<RANK> extents)
{
  static const auto increment=[&](auto n, auto inc)
  {
    using namespace hana::literals;
    if constexpr (n!=0_c) {
      if (idx[n]==extents[n]-1) {idx[n]=0; inc(n-1_c,inc);}
      else idx[n]++;
    }
    else idx[0]++; // This will eventually put the iterator into an illegal state, but this is how every (unchecked) iterator works.
  };
  
  increment(hana::llong_c<RANK-1>,increment);
  
  return idx;
}



/// A non-owning view, should be passed by value.
template <typename T, size_t RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default; MultiArrayView() = default;
  
  MultiArrayView& operator=(const MultiArrayView&) = default; MultiArrayView& operator=(MultiArrayView&&) = default;
  
  /* template <typename... Indices> requires ( sizeof...(Indices)==RANK && ( ... && std::is_convertible_v<Indices,size_t> ) ) T& operator() (Indices... i); */
  
  /// offset can come from a previous slicing
  MultiArrayView(Extents<RANK> e, Extents<RANK> s, size_t o, auto&&... spanArgs)
    : extents{e}, strides{s}, offset{o}, dataView{std::forward<decltype(spanArgs)>(spanArgs)...} {}

  /// implicit conversion to a const view
  operator MultiArrayView<const T, RANK>() requires ( !std::is_const<T>() )
  {
    return {extents,strides,offset};
  }
    
  /// A MultiArray(View) is nothing else than a function that takes a set of indices (=multi-index) as argument and returns the element corresponding to that multi-index
  /** The index of the underlying 1D storage is calculated as 
   * \f[ o + \sum_{\iota=0}^{\iota<R} s_{\iota} i_{\iota} ,\f]
   * where \f$o\f$ is the offset, \f$R\f$ is the rank, \f$s\f$ is the set of strides and \f$i\f$ is the multi-index.
   */
  T& operator() (Extents<RANK> indices) const
  {
#ifndef   NDEBUG
    for (size_t i=0; i<RANK; ++i)
      if (indices[i] >= extents[i])
        throw std::range_error("Index position: "+std::to_string(i)+", index value: "+std::to_string(indices[i])+", extent: "+std::to_string(extents[i]));
#endif // NDEBUG
    return dataView[cppqedutils::ranges::fold(boost::combine(indices,strides),offset,
                                              [&](auto init, auto ids) {return init+ids.template get<0>()*ids.template get<1>();} ) ];
  }
  
  T& operator() (std::convertible_to<size_t> auto ... i) const
  requires (sizeof...(i)==RANK)
  {
    return operator()(Extents<RANK>{i...});
  }
  
  /// A simple specialization for unary views
  T& operator() (std::convertible_to<size_t> auto i) const requires (RANK==1) {return dataView[offset+strides[0]*i];}

  /// A simple specialization for binary views
  T& operator() (std::convertible_to<size_t> auto i, std::convertible_to<size_t> auto j) const requires (RANK==2) {return dataView[offset+strides[0]*i+strides[1]*j];}

  Extents<RANK> extents, strides;

  size_t offset;
  
  std::span<T> dataView; //< the storage can eventually become a kokkos::View
  
};


template <typename T, size_t RANK> requires ( !std::is_const<T>() )
using MultiArrayConstView = MultiArrayView<const T, RANK>;


namespace multiarray {

template <size_t RANK>
auto calculateStrides(Extents<RANK> extents)
{
  Extents<RANK> strides; strides[0]=1;
  for (size_t i=1; i<RANK; ++i) strides[i]=strides[i-1]*extents[i-1];
  return strides;
}

template <size_t RANK>
auto calculateExtent(const Extents<RANK>& extents)
{
  return cppqedutils::ranges::fold(extents,size_t{1},std::multiplies{});
}

} // multiarray


/// Owns data, hence it’s uncopyable (and unassignable), but can be move-constructed (move-assigned) 
template <typename T, size_t RANK>
requires ( !std::is_const<T>() )
class MultiArray : public MultiArrayView<T,RANK>
{
public:
  using StorageType = std::vector<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;
  
  /// initializer is a callable, taking the total size as argument.
  MultiArray(const Extents<RANK>& extents,
             std::function<StorageType(size_t)> initializer=[](size_t s) {
    return StorageType(s); // no braces here, please!!!
  })
    : MultiArrayView<T,RANK>{extents,multiarray::calculateStrides(extents),0},
      data_{initializer(multiarray::calculateExtent(extents))}
  {
    this->dataView=std::span<T>(data_);
  }

  friend void tag_invoke( boost::json::value_from_tag, boost::json::value& jv, const MultiArray<T,RANK>& ma )
  {
    jv = {
      { "extents" , ma.extents },
      //    { "strides" , ma.strides },
      { "data", ma.data_ }
    };
  }

  friend auto tag_invoke( boost::json::value_to_tag< MultiArray<T,RANK> >, const boost::json::value& jv )
  {
    const boost::json::object & obj = jv.as_object();
    return MultiArray<T,RANK> {
      value_to<Extents<RANK>>( obj.at( "extents" ) ),
#ifndef   NDEBUG
      [&] (size_t s) {
        auto ret{value_to<StorageType>( obj.at( "data" ) )};
        if (ret.size() != s) throw std::runtime_error("Mismatch in data size and extents parsed from JSON");
        return ret;        
      }
#else  // NDEBUG
      [&] (size_t) {return value_to<StorageType>( obj.at( "data" ) ); }
#endif // NDEBUG
    };
  }
  
  
private:
  StorageType data_;
  
};


namespace multiarray {

  
/// Checking the consistency of template arguments for use in slicing
/**
 * - size of `RetainedAxes` must not be larger than `RANK`
 * - `RetainedAxes` must not contain duplicated elements
 * - all elements of `RetainedAxes` must be smaller than `RANK`
 * 
 * \todo Think over requirements here & everywhere else in multiarray
 */
constexpr bool consistent(auto retainedAxes)
{/*
  if constexpr ( hana::sequence(retainedAxes) ) {
    constexpr auto sorted = hana::sort(retainedAxes);
    return (hana::unique(sorted) == sorted);    
  }
  else*/ return true;
}


/// Filters in the indices corresponding to a subsystem
/**
 * When the size of `retainedAxes` equals `RANK`, this is a transposition
 */
template <size_t RANK>
auto filterIn(Extents<RANK> idx, auto retainedAxes) requires ( hana::size(retainedAxes) <= RANK )
{
  Extents<hana::size(retainedAxes)> res;
  hana::fold(retainedAxes, res.begin(), [&] (auto iterator, auto ind) {*iterator=idx[ind]; return ++iterator;});
  return res;
}

  
/// Filters out the indices corresponding to a subsystem (specified by `retainedAxes`)
template<size_t RANK>
auto filterOut(Extents<RANK> idx, auto retainedAxes) requires ( hana::size(retainedAxes) <= RANK )
{
  Extents<RANK-hana::size(retainedAxes)> res;
  
  hana::fold(tmptools::ordinals<RANK>, res.begin(), [&] (auto iterator, auto ind) {
    if ( ! hana::contains(retainedAxes, ind) ) *iterator++ = idx[ind];
    return iterator;
  } );
  
  return res;

}


/// Calculates the slices offsets purely from the (dummy) extents and strides
/**
 * Any indexing offset that the MultiArrayView might itself have is added at the point of indexing/slicing
 */
template <size_t RANK>
auto calculateSlicesOffsets(Extents<RANK> extents, Extents<RANK> strides, auto retainedAxes)
{
  Extents<RANK-hana::size(retainedAxes)> 
    dummyExtents{filterOut(extents,retainedAxes)},
    dummyStrides{filterOut(strides,retainedAxes)}, 
    idx; idx.fill(0);
  
  std::vector<size_t> res(calculateExtent(dummyExtents));

  for (auto i=res.begin(); i!=res.end(); (
    *i++ = cppqedutils::ranges::fold( boost::combine(idx,dummyStrides), 
                                      0, 
                                      [&](auto init, auto ids) {return init+ids.template get<0>()*ids.template get<1>();} ),
    incrementMultiIndex(idx,dummyExtents)));
  
  return res;
}


} // multiarray



/// To be used by Composites
template <size_t RANK>
auto calculateSlicesOffsets(Extents<RANK> extents, auto retainedAxes)
{
  return multiarray::calculateSlicesOffsets(extents,calculateStrides(extents),retainedAxes);
}


/**
 * \todo SlicesRange should be a class with deduction guides and a member class Iterator derived from boost::random_access_iterator_helper
 * (Storing full MultiArrayView classes is somewhat wasteful, since extents and strides are the same for all of them in a given SlicesRange)
 */
template <typename T, size_t RANK>
auto slicesRange(MultiArrayView<T,RANK> mav, auto retainedAxes, auto&& offsets)
requires ( RANK>1 && (hana::maximum(retainedAxes) < hana::integral_c<::size_t,RANK>).value && multiarray::consistent(retainedAxes) )
{
  using SliceType=MultiArrayView<T,hana::size(retainedAxes)>;
  
  std::vector<SliceType> res(offsets.size());
  
  for (auto&& [slice,sliceOffset] : boost::combine(res,offsets) )
    slice=SliceType{multiarray::filterIn(mav.extents,retainedAxes),multiarray::filterIn(mav.strides,retainedAxes),mav.offset+sliceOffset,mav.dataView};
  
  return res;
}
  


/// Value equality, so offset doesn’t matter
/** \todo here, the full data is compared, whereas m1 and/or m2 can be only slices */
template <typename T, size_t RANK>
bool operator==(MultiArrayView<T,RANK> m1, MultiArrayView<T,RANK> m2)
{
  return m1.extents==m2.extents && std::ranges::equal(m1.dataView,m2.dataView);
}
 

} // cppqedutils
