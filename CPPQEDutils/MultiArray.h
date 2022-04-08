// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "TMP_Tools.h"

#include <array>
#include <concepts>
#include <functional>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>


namespace cppqedutils {

  
template <unsigned RANK>
using Extents=std::array<size_t,RANK>;


namespace multiarray {

/// Checking the consistency of template arguments for use in slicing
/**
 * - size of `RetainedAxes` must not be larger than `RANK`
 * - `RetainedAxes` must not contain duplicated elements
 * - all elements of `RetainedAxes` must be smaller than `RANK`
 * 
 * \see \ref specifyingsubsystems
 * 
 */
template <int RANK, typename RetainedAxes>
constexpr bool consistent(RetainedAxes ra) // if the optional argument is not given, RetainedAxes needs to be default-constructible 
{
  using namespace hana; using namespace literals;

  constexpr bool res = maximum(ra) < int_c<RANK> && minimum(ra) >= 0_c;
  
  if constexpr ( Sequence<RetainedAxes>::value ) {
    constexpr auto sorted = sort(ra);
    res &= (unique(sorted) == sorted);    
  }
  
  return res;
}

} // multiarray


template <typename T, unsigned RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default;
  
  MultiArrayView& operator=(const MultiArrayView&) = default; MultiArrayView& operator=(MultiArrayView&&) = default;
  
  /* template <typename... Indices> requires ( sizeof...(Indices)==RANK && ( ... && std::is_convertible_v<Indices,size_t> ) ) T& operator() (Indices... i); */
  
  MultiArrayView(const Extents<RANK>& extents, const Extents<RANK>& strides, size_t offset, auto&&... spanArgs)
  : extents_{extents}, strides_{strides}, offset_{offset}, dataView_{std::forward<decltype(spanArgs)>(spanArgs)...} {}
  
  T& operator() (const Extents<RANK>& indices) const
  {
#ifndef   NDEBUG
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {
      if ( size_t actualIndex=indices[idx], actualExtent=extents_[idx]; actualIndex >= actualExtent )
        throw std::range_error("Index position: "+std::to_string(idx)+", index value: "+std::to_string(actualIndex)+", extent: "+std::to_string(actualExtent));
    });
#endif // NDEBUG
    return dataView_[hana::fold(tmptools::ordinals<RANK>,offset_,[&](size_t v,auto idx) {return v+strides_[idx]*indices[idx];})];
    
  }
  
  T& operator() (std::convertible_to<size_t> auto ... i) const
  requires (sizeof...(i)==RANK)
  {
    return operator()(Extents<RANK>{i...});
  }
  
  /// A simple specialization for unary views
  T& operator() (std::convertible_to<size_t> auto i) const requires (RANK==1) {return dataView_[offset_+strides_[0]*i];}

  /// A simple specialization for binary views
  T& operator() (std::convertible_to<size_t> auto i, std::convertible_to<size_t> auto j) const requires (RANK==2) {return dataView_[offset_+strides_[0]*i+strides_[1]*j];}

  
  template <typename RetainedAxes>
  auto slicesRange(RetainedAxes ra, std::vector<size_t>&& offsets) const requires (multiarray::consistent<RANK>(ra))
  {
    std::vector<MultiArrayView<T,hana::size(ra)>> res;
  }
  
 
private:
  Extents<RANK> extents_, strides_;

  size_t offset_;
  
protected:
  std::span<T> dataView_; // the storage can eventually become a kokkos::View
  
};


namespace multiarray {

template <unsigned RANK>
auto calculateStrides(const Extents<RANK>& extents)
{
  Extents<RANK> strides; strides[0]=1;
  for (size_t i=1; i<RANK; ++i) strides[i]=strides[i-1]*extents[i-1];
  return strides;
}

} // multiarray



template <typename T, unsigned RANK>
class MultiArray : public MultiArrayView<T,RANK>
{
public:
  using StorageType = std::vector<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;
  
  template<typename INITIALIZER=std::function<void(StorageType&)>>
  MultiArray(const Extents<RANK>& extents, INITIALIZER&& initializer=[](StorageType&) {})
  : MultiArrayView<T,RANK>{extents,multiarray::calculateStrides<RANK>(extents),0},
  data_(std::accumulate(extents.begin(),extents.end(),1ul,std::multiplies{}))
  {
    this->dataView_=std::span<T>(data_);
    initializer(data_); // initialize vector
  }

  
private:
  StorageType data_;
  
};


namespace multiarray {

template <unsigned RANK, typename RetainedAxes>
auto calculateSlicesOffsets(Extents<RANK> extents, Extents<RANK> strides, RetainedAxes ra)
{
  std::vector<size_t> res;
  return res;
}


template <unsigned RANK, typename RetainedAxes>
auto calculateSlicesOffsets(Extents<RANK> extents, RetainedAxes ra)
{
  return calculateSlicesOffsets(extents,calculateStrides(extents),ra);
}


} // multiarray


} // cppqedutils
