// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Algorithm.h"
#include "TMP_Tools.h"

#include <boost/range/combine.hpp>

#include <array>
#include <concepts>
#include <functional>
#include <span>
#include <stdexcept>
#include <vector>


namespace cppqedutils {

  
template <size_t RANK>
using Extents=std::array<size_t,RANK>;


template <typename T, size_t RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default; MultiArrayView() = default;
  
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

  
  const Extents<RANK>& extents() const {return extents_;}
  const Extents<RANK>& strides() const {return strides_;}
  
  size_t offset() const {return offset_;}
  
  std::span<T> dataView() const {return dataView_;}
  
private:
  Extents<RANK> extents_, strides_;

  size_t offset_;
  
protected:
  std::span<T> dataView_; // the storage can eventually become a kokkos::View
  
};


namespace multiarray {

template <size_t RANK>
auto calculateStrides(const Extents<RANK>& extents)
{
  Extents<RANK> strides; strides[0]=1;
  for (size_t i=1; i<RANK; ++i) strides[i]=strides[i-1]*extents[i-1];
  return strides;
}

template <size_t RANK>
auto calculateExtent(const Extents<RANK>& extents)
{
  return cppqedutils::ranges::fold(extents,1ul,std::multiplies{});
}

} // multiarray



template <typename T, size_t RANK>
class MultiArray : public MultiArrayView<T,RANK>
{
public:
  using StorageType = std::vector<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;
  
  template<typename INITIALIZER=std::function<void(StorageType&)>>
  MultiArray(const Extents<RANK>& extents, INITIALIZER&& initializer=[](StorageType&) {})
  : MultiArrayView<T,RANK>{extents,multiarray::calculateStrides(extents),0}, data_(multiarray::calculateExtent(extents))
  {
    this->dataView_=std::span<T>(data_);
    initializer(data_); // initialize vector
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
 * \see \ref specifyingsubsystems
 * 
 */
template <size_t RANK, typename RetainedAxes>
constexpr bool consistent(RetainedAxes ra) // if the optional argument is not given, RetainedAxes needs to be default-constructible 
{
  using namespace hana;
  if constexpr ( Sequence<RetainedAxes>::value ) {
    constexpr auto sorted = sort(ra);
    return (unique(sorted) == sorted);    
  }
  else return true;
}


template <size_t RANK, typename RetainedAxes>
auto filterIn(const Extents<RANK>& extents, RetainedAxes ra)
{
  Extents<hana::size(ra)> res;
  hana::fold(ra,res.begin(),
             [&] (auto iterator, auto ind) {*iterator=extents[ind]; return ++iterator;});
  return res;
}

  
/// Filters out the indices corresponding to a subsystem
/**
 * \tparamRANK
 * \tparam RetainedAxes compile-time vector specifying the subsystem (cf. \ref specifyingsubsystems)
 * 
 * \param idx the indices to be filtered (of size `RANK`)
 * 
 * \return the indices *not* contained by the subsystem specified by `RetainedAxes` (dummy indices)
 * 
 */
template<size_t RANK, typename RetainedAxes>
auto filterOut(const Extents<RANK>& idx, RetainedAxes ra)
{
  using namespace hana;
  
  ::size_t origIndex=0, resIndex=0;
  
  static constexpr auto sorted = [&]() {
    if constexpr ( Sequence<RetainedAxes>::value ) return sort(ra);  
    else return ra;
  }();
  
  Extents<RANK-hana::size(ra)> res;

  for_each(sorted,[&](auto i) {
    for (; origIndex<i; (++origIndex,++resIndex)) res[resIndex]=idx[origIndex];
    ++origIndex; // skip value found in sorted
  });
  // the final segment:
  for (; origIndex<idx.size(); (++origIndex,++resIndex)) res[resIndex]=idx[origIndex];
  
  return res;

}


/// This calculates the slices offsets purely from the (dummy) extents and strides
/**
 * Any indexing offset that the MultiArrayView might itself have is added at the point of indexing/slicing
 */
template <size_t RANK, typename RetainedAxes>
auto calculateSlicesOffsets(Extents<RANK> extents, Extents<RANK> strides, RetainedAxes ra)
{
  Extents<RANK-hana::size(ra)> dummyExtents{filterOut(extents,ra)}, dummyStrides{filterOut(strides,ra)}, 
                               idx; idx.fill(0);
  
  std::vector<size_t> res(calculateExtent(dummyExtents));

  auto increment=[&](auto n, auto inc)
  {
    using namespace hana; using namespace literals;
    if constexpr (n!=0_c) {
      if (idx[n]==dummyExtents[n]-1) {idx[n]=0; inc(n-1_c,inc);}
      else idx[n]++;
    }
    else idx[0]++; // This will of course eventually put the iterator into an illegal state when idx(0)>ubound(0), but this is how every (unchecked) iterator works.
  };
  
  for (auto i=res.begin(); i!=res.end(); (
    *i++ = cppqedutils::ranges::fold( boost::combine(idx,dummyStrides), 
                                      0, 
                                      [&](auto init, auto ids) {return init+ids.template get<0>()*ids.template get<1>();} ),
    increment(hana::llong_c<idx.size()-1>,increment)));
  
  return res;
}


} // multiarray



/// To be used by Composites
template <size_t RANK, typename RetainedAxes>
auto calculateSlicesOffsets(Extents<RANK> extents, RetainedAxes ra)
{
  return multiarray::calculateSlicesOffsets(extents,calculateStrides(extents),ra);
}


/// TODO: SlicesRange should be a class with deduction guides and a member class Iterator derived from boost::random_access_iterator_helper
template <typename T, size_t RANK, typename RetainedAxes>
auto slicesRange(const MultiArrayView<T,RANK>& mav, RetainedAxes ra, std::vector<size_t>&& offsets)
requires ( RANK>1 && (hana::maximum(ra) < hana::integral_c<::size_t,RANK>).value && multiarray::consistent<RANK>(ra) )
{
  using SliceType=MultiArrayView<T,hana::size(ra)>;
  
  std::vector<SliceType> res(offsets.size());
  
  for (auto&& [slice,sliceOffset] : boost::combine(res,offsets) )
    slice=SliceType{multiarray::filterIn(mav.extents(),ra),multiarray::filterIn(mav.strides(),ra),mav.offset()+sliceOffset,mav.dataView()};
  
  return res;
}
  


} // cppqedutils
