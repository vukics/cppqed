// Copyright András Vukics 2022–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Algorithm.h"
#include "Archive.h"
#include "Traits.h"

#include <boost/serialization/vector.hpp>

#include <array>
#include <concepts>
#include <span>
#include <stdexcept>
//#include <valarray>
#include <vector>


namespace cppqedutils {

/// should be passed by value, but can be `std::move`d when used in class constructors
template <size_t RANK>
using Extents=std::array<size_t,RANK>;


/// endpoint guard is: `idx[0]!=extents[0]` note the `[0]`!
/// TODO: this could be formulated as a std compliant range, cf. chatGPT’s implementation towards the end of SliceIterator.h
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
    else idx[0]++; // This will eventually put the index into an illegal state, but this is how all (unchecked) iterators work.
  };
  
  increment(hana::llong_c<RANK-1>,increment);
  
  return idx;
}

inline auto& incrementMultiIndex(Extents<0>& idx, Extents<0>) {return idx;}


namespace multiarray {

static constexpr struct Placeholder {} _;

} // multiarray


/// A non-owning view, should be passed by value.
/**
 * Follows reference semantics regarding constness.
 * All four possibilities of (const)MultiArrayView<(const)T,RANK> can be used
 * and have the same semantics as if MultiArrayView was a pointer (same as with std::span)
 *
 * TODO: alternative design: MultiArrayBase with StorageType as template parameter
 * then, when StorageType is a span, it’s a view, and when it’s a vector, it’s an owning MultiArray.
 * In this case, MultiArray need not be derived from MultiArrayView, of course :)
 * In this case, StorageType could even be a std range view (e.g. a transformed view),
 * in which case MultiArrayView could represent the result of lazy operations
 */
template <typename T, size_t RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default; MultiArrayView() = default;
  MultiArrayView& operator=(const MultiArrayView&) = default; MultiArrayView& operator=(MultiArrayView&&) = default;
  
  /// offset can come from a previous slicing
  MultiArrayView(Extents<RANK> extents, Extents<RANK> strides, size_t offset, auto&&... dataView)
    : extents{extents}, strides{strides}, offset{offset}, dataView{std::forward<decltype(dataView)>(dataView)...} {}

  /// implicit conversion to a const view
  operator MultiArrayView<const T, RANK>() const requires ( !std::is_const_v<T> ) {return {extents,strides,offset,dataView}; }

  void checkBounds(std::convertible_to<size_t> auto ... i) const
  {
#ifndef   NDEBUG
    auto e=extents.begin();
    (... , [&] {
      if (auto eComp=e++; i >= *eComp)
        throw std::range_error("Index position: "+std::to_string(eComp-extents.begin())+", index value: "+std::to_string(i)+", extent: "+std::to_string(*eComp));
    } () );
#endif // NDEBUG
  }

  /// A MultiArray(View) is nothing else than a function that takes a set of indices (=multi-index) as argument and returns the element corresponding to that multi-index
  /** 
   * The index of the underlying 1D storage is calculated as \f[ o + \sum_{\iota=0}^{\iota<R} s_{\iota} i_{\iota} ,\f]
   * where \f$o\f$ is the offset, \f$R\f$ is the rank, \f$s\f$ is the set of strides and \f$i\f$ is the multi-index.
   * 
   * Note that the subscription operator is always const, as it doesn’t modify the view itself, only maybe the underlying data.
   */
  T& operator() (Extents<RANK> idx) const
  {
    std::apply([this] (auto &&... args) { checkBounds(std::forward<decltype(args)>(args)...); },idx);
    return dataView[std::ranges::fold_left( std::views::zip(idx,strides), offset,
                                            [&] (auto init, auto ids) {return init+get<0>(ids)*get<1>(ids);} ) ];
  }
  
  T& operator() (std::convertible_to<size_t> auto ... i) const requires (sizeof...(i)==RANK)
  {
    checkBounds(i...);
    size_t idx=0;
    auto s=strides.begin();
    return dataView[ offset + ( ... , [&] { return idx+=(*s++)*i;} () ) ];
  }
  
  /// A simple specialization for unary views
  T& operator() (std::convertible_to<size_t> auto i) const requires (RANK==1) {checkBounds(i); return dataView[offset+strides[0]*i];}

  /// A simple specialization for binary views
  T& operator() (std::convertible_to<size_t> auto i, std::convertible_to<size_t> auto j) const requires (RANK==2) {
    checkBounds(i,j); return dataView[offset+strides[0]*i+strides[1]*j];
  }

  /// TODO: extend these to arbitrary RANK
  auto operator() (std::convertible_to<size_t> auto i, multiarray::Placeholder) const requires (RANK==2) {
    return MultiArrayView<T,1>{ {extents[1]}, {strides[1]}, offset+i*strides[0], dataView};
  }

  auto operator() (multiarray::Placeholder, std::convertible_to<size_t> auto j) const requires (RANK==2) {
    return MultiArrayView<T,1>{ {extents[0]}, {strides[0]}, offset+j*strides[1], dataView};
  }

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


template <typename T1, typename T2, size_t RANK>
void checkExtents(MultiArrayConstView<T1,RANK> m1, MultiArrayConstView<T2,RANK> m2, std::string message)
{
#ifndef   NDEBUG
  if (m1.extents!=m2.extents) throw std::runtime_error("Extent mismatch in "+message+": "+json(m1.extents).dump()+" "+json(m2.extents).dump());
#endif // NDEBUG
}

/// Element-by-element comparison
/** \note The operator-style syntax is not used, since this can be a rather expensive operation! */
template <typename T1, typename T2, size_t RANK>
bool isEqual(MultiArrayConstView<T1,RANK> m1, MultiArrayConstView<T2,RANK> m2)
{
  checkExtents(m1,m2,"MultiArrayView comparison");
  bool res=true;
  for (Extents<RANK> idx{}; idx[0]!=m1.extents[0]; incrementMultiIndex(idx,m1.extents))
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
  return std::ranges::fold_left(extents,1uz,std::multiplies{});
}

template <typename T>
const auto noInit = [] (size_t e) -> std::vector<T> {return std::vector<T>(e);};

template <typename T>
const auto zeroInit = [] (size_t e) -> std::vector<T> {return std::vector<T>(e,T(0.));};


template <typename T>
using StorageType = std::vector<T>;


} // multiarray


/// Owns data, hence it’s uncopyable (and unassignable), but can be move-constructed (move-assigned)
/** Follows value semantics regarding constness. */
template <typename T, size_t RANK>
requires ( !std::is_const_v<T> )
class MultiArray : public MultiArrayConstView<T,RANK>
{
public:
  /// the storage might eventually become a kokkos::View or a std::valarray
  /** TODO: check whether `std::vector` incurs too much overhead – with std::valarray, the slicing functionality could be exploited
   * (MultiArrayView could store a std::slice_array in this case? Ooops, std::slice_array cannot be further sliced???!!! :-O )
   */
  using StorageType = multiarray::StorageType<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;

  /// initializer is a callable, taking the total size as argument.
  MultiArray(Extents<RANK> extents, auto&& initializer) requires requires (size_t e) { { initializer(e) } -> std::convertible_to<StorageType>;}
    : MultiArrayConstView<T,RANK>{extents,multiarray::calculateStrides(extents),0}, data_{initializer(multiarray::calculateExtent(extents))}
  {
    this->dataView=std::span<T>(data_);
  }

  explicit MultiArray(Extents<RANK> extents) : MultiArray{extents,multiarray::noInit<T>} {}

  friend MultiArray copy(const MultiArray& ma) {return MultiArray{ma.extents, [&] (size_t) {return ma.dataStorage();}};}

  /// Element-by-element assignment
  /** \note The operator-style syntax is not used, since this can be a rather expensive operation! */
  template <typename TT>
  void assignTo(MultiArrayConstView<TT,RANK> macv)
  {
#ifndef   NDEBUG
    if (this->dataView.data()==macv.dataView.data())
      throw std::runtime_error("Self assignment attempted in MultiArray");
#endif // NDEBUG
    checkExtents(*this,macv,"MultiArrayView assignment");

    for (Extents<RANK> idx{}; idx[0]!=this->extents[0]; incrementMultiIndex(idx,this->extents)) (*this)(idx)=macv(idx);
  }

  /// in contexts where implicit conversions are ignored, it is sometimes convenient to have this as a member :)
  MultiArrayConstView<T,RANK> constView() const {return *this;}
  
  /// conversion to a view
  auto mutableView() {return MultiArrayView<T, RANK>{this->extents,this->strides,0,data_};}
  
  operator MultiArrayView<T,RANK>() {return mutableView();}

  /// non-const subscripting
  T& operator()(auto&&... i) {return const_cast<T&>(static_cast<MultiArrayConstView<T,RANK>>(*this)(std::forward<decltype(i)>(i)...)) ;}
  
  /// this overload is needed, otherwise any call to operator() with a MultiArray object resolves to the previous function
  const T& operator()(auto&&... i) const {return static_cast<MultiArrayConstView<T,RANK>>(*this)(std::forward<decltype(i)>(i)...) ;}

  /// access data storage
  StorageType& dataStorage() {return data_;}
  const StorageType& dataStorage() const {return data_;}

  /// JSONize
  friend void to_json( json& jv, const MultiArray& ma )
  {
    jv = json{
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
  
private:
  StorageType data_;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive& ar, const unsigned int) const {ar & this->extents & data_;}

  template<class Archive> void load(Archive& ar, const unsigned int)
  {
    ar & this->extents & data_; this->strides=multiarray::calculateStrides(this->extents); this->offset=0;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
};

template <typename T, size_t RANK>
constexpr auto passByValue_v<MultiArray<T,RANK>> = false;


namespace multiarray {


template <typename T, size_t RANK>
const auto copyInit(const std::ranges::sized_range auto& input) {return [&] (size_t e) {
  if (size(input) != e) throw std::runtime_error("Extent mismatch in MultiArray copyInit: "+std::to_string(size(input))+" "+std::to_string(e));
  return typename MultiArray<T,RANK>::StorageType(begin(input),end(input));
};}


} // multiarray


} // cppqedutils
