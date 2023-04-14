// Copyright András Vukics 2022–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Algorithm.h"
#include "Archive.h"
#include "Traits.h"

#include <boost/operators.hpp>

#include <boost/range/combine.hpp>

#include <array>
#include <concepts>
#include <span>
#include <stdexcept>
//#include <valarray>
#include <vector>


namespace cppqedutils {


template<size_t... i>
constexpr auto retainedAxes = hana::tuple_c<size_t,i...>;

template<size_t begin, size_t end>
constexpr auto compileTimeRange = hana::range_c<size_t,begin,end>;

template<size_t end>
constexpr auto compileTimeOrdinals = compileTimeRange<0,end>;

/// should be passed by value
template <size_t RANK>
using Extents=std::array<size_t,RANK>;


/// TODO: this could be formulated as a std compliant range with `extents` as an endpoint guard
template <size_t RANK>
auto& incrementMultiIndex(Extents<RANK>& idx, Extents<RANK> extents)
{
  const auto increment=[&](auto n, auto inc)
  {
    using namespace hana::literals;
    if constexpr (n!=0_c) {
      if (idx[n]==extents[n]-1) {idx[n]=0; inc(n-1_c,inc);}
      else idx[n]++;
    }
    else idx[0]++; // This will eventually put the iterator into an illegal state, but this is how all (unchecked) iterators work.
  };
  
  increment(hana::llong_c<RANK-1>,increment);
  
  return idx;
}



/// A non-owning view, should be passed by value.
/**
 * Follows reference semantics regarding constness.
 * All four possibilities of (const)MultiArrayView<(const)T,RANK> can be used and have the same semantics as if MultiArrayView was a pointer (same as with std::span)
 */
template <typename T, size_t RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default; MultiArrayView() = default;
  MultiArrayView& operator=(const MultiArrayView&) = default; MultiArrayView& operator=(MultiArrayView&&) = default;
  
  /* template <typename... Indices> requires ( sizeof...(Indices)==RANK && ( ... && std::is_convertible_v<Indices,size_t> ) ) T& operator() (Indices... i); */
  
  /// offset can come from a previous slicing
  MultiArrayView(Extents<RANK> e, Extents<RANK> s, size_t o, auto&&... dv)
    : extents{e}, strides{s}, offset{o}, dataView{std::forward<decltype(dv)>(dv)...} {}

  /// implicit conversion to a const view
  operator MultiArrayView<const T, RANK>() const requires ( !std::is_const_v<T> ) {return {extents,strides,offset,dataView}; }

  /// A MultiArray(View) is nothing else than a function that takes a set of indices (=multi-index) as argument and returns the element corresponding to that multi-index
  /** 
   * The index of the underlying 1D storage is calculated as \f[ o + \sum_{\iota=0}^{\iota<R} s_{\iota} i_{\iota} ,\f]
   * where \f$o\f$ is the offset, \f$R\f$ is the rank, \f$s\f$ is the set of strides and \f$i\f$ is the multi-index.
   * 
   * Note that the subscription operator is always const, as it doesn’t modify the view itself, only maybe the underlying data.
   */
  T& operator() (Extents<RANK> idx) const
  {
#ifndef   NDEBUG
    for (size_t i=0; i<RANK; ++i) if (idx[i] >= extents[i]) throw std::range_error("Index position: "+std::to_string(i)+", index value: "+std::to_string(idx[i])+", extent: "+std::to_string(extents[i]));
#endif // NDEBUG
    return dataView[std_ext::ranges::fold(boost::combine(idx,strides),offset,
                                          [&](auto init, auto ids) {return init+ids.template get<0>()*ids.template get<1>();} ) ];
  }
  
  
  T& operator() (std::convertible_to<size_t> auto ... i) const requires (sizeof...(i)==RANK)
  {
    size_t idx=0;
    auto e=extents.begin(), s=strides.begin();

    return dataView[ offset + ( ... , [&] {
#ifndef   NDEBUG
      if (auto eComp=e++; i >= *eComp) throw std::range_error("Index position: "+std::to_string(eComp-extents.begin())+", index value: "+std::to_string(i)+", extent: "+std::to_string(*eComp));
#endif // NDEBUG
      return idx+=(*s++)*i;} () ) ];

  }
  
  /// A simple specialization for unary views
  T& operator() (std::convertible_to<size_t> auto i) const requires (RANK==1) {return dataView[offset+strides[0]*i];}

  /// A simple specialization for binary views
  T& operator() (std::convertible_to<size_t> auto i, std::convertible_to<size_t> auto j) const requires (RANK==2) {return dataView[offset+strides[0]*i+strides[1]*j];}

  Extents<RANK> extents, strides;

  size_t offset;
  
  std::span<T> dataView;
  
};


template <typename T, size_t RANK>
constexpr auto passByValue_v<MultiArrayView<T,RANK>> = true;


template <typename T, size_t RANK> requires ( !std::is_const_v<T> )
using MultiArrayConstView = MultiArrayView<const T, RANK>;


// MultiArrayView is itself a reference
template <typename T, size_t RANK> requires ( !std::is_const_v<T> )
struct ReferenceMF<MultiArrayView<T,RANK>> : std::type_identity<MultiArrayView<T,RANK>> {};

template <typename T, size_t RANK> requires ( !std::is_const_v<T> )
struct ConstReferenceMF<MultiArrayView<T,RANK>> : std::type_identity<MultiArrayConstView<T,RANK>> {};



/// Element-by-element comparison
/** \note The operator-style syntax is not used, since this can be a rather expensive operation! */
template <typename T1, typename T2, size_t RANK>
bool isEqual(MultiArrayConstView<T1,RANK> m1, MultiArrayConstView<T2,RANK> m2)
{
#ifndef   NDEBUG
  if (m1.extents!=m2.extents) throw std::runtime_error("Extent mismatch in MultiArrayView comparison: "+json(m1.extents).dump()+" "+json(m2.extents).dump());
#endif // NDEBUG
  bool res=true;
  for (Extents<RANK> idx{}; idx!=m1.extents; incrementMultiIndex(idx,m1.extents))
    res &= ( m1(idx)==m2(idx) );
  return res;
}


namespace multiarray {

template <size_t RANK>
auto calculateStrides(Extents<RANK> extents)
{
  Extents<RANK> strides; strides[0]=1;
  for (size_t i=1; i<RANK; ++i) strides[i]=strides[i-1]*extents[i-1];
  return strides;
}

template <size_t RANK>
auto calculateExtent(Extents<RANK> extents)
{
  return std_ext::ranges::fold(extents,size_t{1},std::multiplies{});
}

} // multiarray


/// Owns data, hence it’s uncopyable (and unassignable), but can be move-constructed (move-assigned) 
/** Follows value semantics regarding constness. */
template <typename T, size_t RANK, template <typename, size_t> typename BASE=MultiArrayConstView>
requires ( !std::is_const_v<T> )
class MultiArray : public BASE<T,RANK>
{
public:
  /// the storage might eventually become a kokkos::View or a std::valarray
  /** TODO: check whether `std::vector` incurs too much overhead – with std::valarray, the slicing functionality could be exploited
   * (MultiArrayView could store a std::slice_array in this case? Ooops, std::slice_array cannot be further sliced???!!! :-O )
   */
  using StorageType = std::vector<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;
  
  /// initializer is a callable, taking the total size as argument.
  MultiArray(Extents<RANK> extents, auto&& initializer) requires requires (Extents<RANK> e) {initializer(e);}
    : BASE<T,RANK>{extents,multiarray::calculateStrides(extents),0}, data_{initializer(extents)}
  {
    this->dataView=std::span<T>(data_);
  }

  explicit MultiArray(Extents<RANK> extents) : MultiArray{extents,[](Extents<RANK> extents) {
    return StorageType(multiarray::calculateExtent(extents)); /* no braces here, please!!! */
  }} {}

  /// Element-by-element assignment
  /** \note The operator-style syntax is not used, since this can be a rather expensive operation! */
  template <typename TT>
  void assignTo(MultiArrayConstView<TT,RANK> macv)
  {
#ifndef   NDEBUG
    if (this->dataView.data()==macv.dataView.data())
      throw std::runtime_error("Self assignment attempted in MultiArray");
    if (this->extents!=macv.extents)
      throw std::runtime_error("Extent mismatch in MultiArrayView assignment: "+toStringJSON(this->extents)+" "+toStringJSON(macv.extents));
#endif // NDEBUG
    for (Extents<RANK> idx{}; idx!=this->extents; incrementMultiIndex(idx,this->extents)) (*this)(idx)=macv(idx);
  }

  /// conversion to a view
  auto mutableView() {return MultiArrayView<T, RANK>{this->extents,this->strides,0,data_};}
  operator MultiArrayView<T,RANK>() {return mutableView();}

  /// non-const subscripting
  T& operator()(auto&&... i) {return const_cast<T&>(static_cast<BASE<T,RANK>>(*this)(std::forward<decltype(i)>(i)...)) ;}
  
  /// JSONize
  friend void to_json( json& jv, const MultiArray& ma )
  {
    jv = {
      { "extents" , ma.extents },
      //    { "strides" , ma.strides },
      { "data", ma.data_ }
    };
  }

  /// unJSONize
  friend void from_json( const json& jv, MultiArray& ma )
  {
    ma.extents=jv["extents"];
    StorageType dataTemp=jv["data"];
#ifndef   NDEBUG
    if (ma.data_.size() != multiarray::calculateExtent(ma.extents)) throw std::runtime_error("Mismatch in data size and extents parsed from JSON");
#endif // NDEBUG
    ma.data_.swap(dataTemp); ma.dataView=std::span<T>(ma.data_);
    ma.strides=multiarray::calculateStrides(ma.extents);
  }
  
protected:
  StorageType data_;
  
private:
  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive& ar, const unsigned int) const {ar & this->extents & data_;}

  template<class Archive> void load(Archive& ar, const unsigned int) {ar & this->extents & data_; this->strides=calculateStrides(this->extents); this->offset=0;}

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
};


template <typename T, size_t RANK>
constexpr auto passByValue_v<MultiArray<T,RANK>> = false;


namespace multiarray {

  
/// Checking the consistency of template arguments for use in slicing
/**
 * - size of `RetainedAxes` must not be larger than `RANK`
 * - `RetainedAxes` must not contain duplicated elements
 * - all elements of `RetainedAxes` must be smaller than `RANK`
 * 
 * \todo Think over requirements here & everywhere else in multiarray
 */
template <auto retainedAxes>
constexpr bool consistentBase = 
/*{
  if constexpr ( hana::sequence(retainedAxes) ) {
    constexpr auto sorted = hana::sort(retainedAxes);
    return (hana::unique(sorted) == sorted);    
  }
  else*/ true;
//}

template <auto retainedAxes, size_t RANK>
constexpr bool consistent = ( RANK>1 && (hana::maximum(retainedAxes) < hana::integral_c<::size_t,RANK>).value && consistentBase<retainedAxes> );


/// Filters in the indices corresponding to a subsystem
/**
 * When the size of `retainedAxes` equals `RANK`, this is a transposition
 * 
 * \todo accept several arrays at the same time (since this is how it is used by clients)
 */
template <auto retainedAxes, size_t RANK>
auto filterIn(Extents<RANK> idx) requires ( hana::size(retainedAxes) <= RANK )
{
  Extents<hana::size(retainedAxes)> res;
  hana::fold(retainedAxes, res.begin(), [&] (auto iterator, auto ind) {*iterator=idx[ind]; return ++iterator;});
  return res;
}

  
/// Filters out the indices corresponding to a subsystem (specified by `retainedAxes`)
/** \todo accept several arrays at the same time (since this is how it is used by clients) */
template <auto retainedAxes, size_t RANK>
auto filterOut(Extents<RANK> idx) requires ( hana::size(retainedAxes) <= RANK )
{
  Extents<RANK-hana::size(retainedAxes)> res;
  
  hana::fold(compileTimeOrdinals<RANK>, res.begin(), [&] (auto iterator, auto ind) {
    if ( ! hana::contains(retainedAxes, ind) ) *iterator++ = idx[ind];
    return iterator;
  } );
  
  return res;

}


/// Calculates the slices offsets purely from the (dummy) extents and strides
/**
 * Any indexing offset that the MultiArrayView might itself have is added at the point of indexing/slicing
 */
template <auto retainedAxes, size_t RANK>
auto calculateSlicesOffsets(Extents<RANK> extents, Extents<RANK> strides)
{
  Extents<RANK-hana::size(retainedAxes)> 
    dummyExtents{filterOut<retainedAxes>(extents)},
    dummyStrides{filterOut<retainedAxes>(strides)}, 
    idx; idx.fill(0);
  
  std::vector<size_t> res(calculateExtent(dummyExtents));

  for (auto i=res.begin(); i!=res.end(); (
    *i++ = std_ext::ranges::fold( boost::combine(idx,dummyStrides), 
                                  0, 
                                  [](size_t init, auto ids) {return init+ids.template get<0>()*ids.template get<1>();} ),
    incrementMultiIndex(idx,dummyExtents)));
    // note: the appearance of idx on both sides of the comma operator cannot cause problems since it ensures left-to-right evaluation order

  return res;
}


} // multiarray



/// To be used by Composites
template <auto retainedAxes, size_t RANK>
auto calculateSlicesOffsets(Extents<RANK> extents)
{
  return multiarray::calculateSlicesOffsets<retainedAxes>(extents,multiarray::calculateStrides(extents));
}



template <auto retainedAxes, typename T, size_t RANK>
auto sliceRangeSimple(MultiArrayView<T,RANK> mav, const std::vector<size_t>& offsets)
requires multiarray::consistent<retainedAxes,RANK>
{
  using SliceType=MultiArrayView<T,hana::size(retainedAxes)>;
  
  std::vector<SliceType> res(offsets.size());
  
  for (auto&& [slice,sliceOffset] : boost::combine(res,offsets) )
    slice=SliceType{multiarray::filterIn<retainedAxes>(mav.extents),multiarray::filterIn<retainedAxes>(mav.strides),mav.offset+sliceOffset,mav.dataView};
  
  return res;
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRangeSimple(MultiArrayView<T,RANK> mav)
{
  return sliceRangeSimple<retainedAxes>(mav,multiarray::calculateSlicesOffsets<retainedAxes>(mav.extents,mav.strides));
}



template <typename Offsets, typename T, size_t RRANK> // RRANK stands for retained rank
requires (std::is_same_v<Offsets,std::span<const size_t>> || std::is_same_v<Offsets,std::vector<size_t>>)
class SliceIterator : public boost::forward_iterator_helper<SliceIterator<Offsets,T,RRANK>,MultiArrayView<T,RRANK>>
/// \todo Unfortunately, std::iterator has been deprecated as of C++17, and there doesn’t seem to be a replacement in STL
{
public:
  using member_iterator = std::conditional_t<std::is_same_v<Offsets,std::span<const size_t>>,std::span<const size_t>::iterator,std::vector<size_t>::const_iterator>;
  
  SliceIterator(Extents<RRANK> e, Extents<RRANK> s, size_t originalOffset, member_iterator i, std::span<T> dv) :
    iter{i}, originalOffset_{originalOffset}, mav_{e,s,*i,dv} {}
  
  auto& operator*() const {return mav_;}
  
  SliceIterator& operator++() {mav_.offset=*(++iter)+originalOffset_; return *this;}
  
  member_iterator iter;

  friend bool operator==(SliceIterator i0, member_iterator i1) {return i0.iter==i1;}
  
private:
  const size_t originalOffset_;
  
  mutable MultiArrayView<T,RRANK> mav_; // the reason why we store a full MultiArrayView is simply that operator* can return a reference
  
};



template <typename Offsets, typename T, size_t RRANK> // RRANK stands for retained rank
class SliceRangeBase
{
public:
  using iterator=SliceIterator<Offsets,T,RRANK>;
  
/*  SliceRangeBase(Extents<RRANK> e, Extents<RRANK> s, Offsets os, std::span<T> dv) requires std::is_same_v<Offsets,std::span<const size_t>>
    : begin{e,s,os.begin(),dv}, end{os.end()}, offsets_{os} {}*/

  SliceRangeBase(Extents<RRANK> e, Extents<RRANK> s, size_t originalOffset, const Offsets& os, std::span<T> dv) // requires std::is_same_v<Offsets,std::vector<size_t>>
    : offsets_{os}, b_{e,s,originalOffset,offsets_.begin(),dv}, e_{offsets_.end()} {}
    
  auto begin() const {return b_;}
  auto end() const {return e_;}
    
private:
  const Offsets offsets_;

  iterator b_;
  
  typename iterator::member_iterator e_;
  
};


/**
 * For the case when Offsets are stored outside of the class.
 * E.g. Composite might calculate and store the offsets for all its subsystems
 */
template <typename T, size_t RRANK>
using SliceRangeReferencing = SliceRangeBase<std::span<const size_t>,T,RRANK>;


/** For the case when Offsets are stored within the class. */
template <typename T, size_t RRANK> // RRANK stands for retained rank
using SliceRangeOwning = SliceRangeBase<std::vector<size_t>,T,RRANK>;



template <auto retainedAxes, typename T, size_t RANK>
auto sliceRange(MultiArrayView<T,RANK> mav, const std::vector<size_t>& offsets) requires multiarray::consistent<retainedAxes,RANK>
{
  return SliceRangeReferencing<T,hana::size(retainedAxes)>{
    multiarray::filterIn<retainedAxes>(mav.extents),
    multiarray::filterIn<retainedAxes>(mav.strides),
    mav.offset,
    std::span<const size_t>(offsets.begin(),offsets.end()),
    mav.dataView};
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRange(MultiArrayView<T,RANK> mav, std::vector<size_t>&& offsets) requires multiarray::consistent<retainedAxes,RANK>
{
  return SliceRangeOwning<T,hana::size(retainedAxes)>{
    multiarray::filterIn<retainedAxes>(mav.extents),
    multiarray::filterIn<retainedAxes>(mav.strides),
    mav.offset,
    std::move(offsets),
    mav.dataView};
}


template <auto retainedAxes, typename T, size_t RANK>
auto sliceRange(MultiArrayView<T,RANK> mav) requires multiarray::consistent<retainedAxes,RANK>
{
  return sliceRange<retainedAxes>(mav,multiarray::calculateSlicesOffsets<retainedAxes>(mav.extents,mav.strides));
}


template <auto axes, typename T, size_t RANK>
MultiArrayView<T,RANK> transpose(MultiArrayView<T,RANK> mav) requires ( hana::size(axes) == RANK )
{
  mav.extents=filterIn<axes>(mav.extents); mav.strides=filterIn<axes>(mav.strides);
  return mav;
}


} // cppqedutils
