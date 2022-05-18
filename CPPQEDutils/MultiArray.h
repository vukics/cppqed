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


namespace cppqedutils {

  
template <size_t RANK>
using Extents=std::array<size_t,RANK>;


/// A non-owning view, should be passed by value.
template <typename T, size_t RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default; MultiArrayView() = default;
  
  MultiArrayView& operator=(const MultiArrayView&) = default; MultiArrayView& operator=(MultiArrayView&&) = default;
  
  /* template <typename... Indices> requires ( sizeof...(Indices)==RANK && ( ... && std::is_convertible_v<Indices,size_t> ) ) T& operator() (Indices... i); */
  
  MultiArrayView(const Extents<RANK>& e, const Extents<RANK>& s, size_t o, auto&&... spanArgs)
  : extents{e}, strides{s}, offset{o}, dataView{std::forward<decltype(spanArgs)>(spanArgs)...} {}
  
  T& operator() (const Extents<RANK>& indices) const
  {
#ifndef   NDEBUG
    for (size_t i=0; i<RANK; ++i)
      if ( size_t actualIndex=indices[i], actualExtent=extents[i]; actualIndex >= actualExtent )
        throw std::range_error("Index position: "+std::to_string(i)+", index value: "+std::to_string(actualIndex)+", extent: "+std::to_string(actualExtent));
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
  
  std::span<T> dataView; // the storage can eventually become a kokkos::View
  
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


/// Owns data, hence it’s uncopyable (and unassignable), but can be move-constructed (move-assigned) 
template <typename T, size_t RANK>
class MultiArray : public MultiArrayView<T,RANK>
{
public:
  using StorageType = std::vector<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;
  
  /// INITIALIZER is a callable, taking Extents<RANK> as argument.
  template<typename INITIALIZER=std::function<StorageType(size_t)>>
  MultiArray(const Extents<RANK>& extents, INITIALIZER&& initializer=[](size_t s) {
    return StorageType(s);
  })
  : MultiArrayView<T,RANK>{extents,multiarray::calculateStrides(extents),0}, data_{initializer(multiarray::calculateExtent(extents))}
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
 */
template <size_t RANK, typename RetainedAxes>
constexpr bool consistent(RetainedAxes ra) // if the optional argument is not given, RetainedAxes needs to be default-constructible 
{
  if constexpr ( hana::Sequence<RetainedAxes>::value ) {
    constexpr auto sorted = hana::sort(ra);
    return (hana::unique(sorted) == sorted);    
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
    using namespace hana::literals;
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


/**
 * TODO: SlicesRange should be a class with deduction guides and a member class Iterator derived from boost::random_access_iterator_helper
 * (Storing full MultiArrayView classes is somewhat wasteful, since extents and strides is the same for all of them in a given SlicesRange)
 */
template <typename T, size_t RANK, typename RetainedAxes, typename OFFSETS>
auto slicesRange(MultiArrayView<T,RANK> mav, RetainedAxes ra, OFFSETS&& offsets)
requires ( RANK>1 && (hana::maximum(ra) < hana::integral_c<::size_t,RANK>).value && multiarray::consistent<RANK>(ra) )
{
  using SliceType=MultiArrayView<T,hana::size(ra)>;
  
  std::vector<SliceType> res(offsets.size());
  
  for (auto&& [slice,sliceOffset] : boost::combine(res,offsets) )
    slice=SliceType{multiarray::filterIn(mav.extents,ra),multiarray::filterIn(mav.strides,ra),mav.offset+sliceOffset,mav.dataView};
  
  return res;
}
  


/// Value equality, so offset doesn’t matter
template <typename T, size_t RANK>
bool operator==(MultiArrayView<T,RANK> m1, MultiArrayView<T,RANK> m2)
{
  return m1.extents==m2.extents && std::ranges::equal(m1.dataView,m2.dataView);
}
 

} // cppqedutils
