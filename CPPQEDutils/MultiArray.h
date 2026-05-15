// Copyright András Vukics 2022–2026. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "Traits.h"

#include <boost/serialization/vector.hpp>

#include <array>
#include <concepts>
#include <span>
#include <stdexcept>
#include <ranges>
#include <vector>


namespace cppqedutils {

/**
 * @brief Concatenates multiple `std::array` instances of different sizes into one.
 *
 * Overcomes the limitation that `std::views::join` cannot join `std::array` instances
 * of different sizes, since the size is part of the type. Uses C++17 fold expressions.
 *
 * @tparam Type The element type, common to all input arrays.
 * @tparam sizes The sizes of the input arrays, deduced from the arguments.
 * @param arrays The arrays to concatenate.
 * @return A new `std::array<Type, (sizes+...)>` containing all elements in order.
 *
 * @see http://stackoverflow.com/a/42774523/1171157
 */
template <typename Type, std::size_t... sizes>
constexpr auto concatenate(const std::array<Type, sizes>&... arrays)
{
  std::array<Type, (sizes + ...)> result;
  std::size_t index{};

  ((std::copy_n(arrays.begin(), sizes, result.begin() + index), index += sizes), ...);

  return result;
}


/// @brief Compile-time–fixed shape and stride descriptor for a rank-`RANK` array. Should be passed by value, but can be `std::move`d in class constructors.
template <size_t RANK>
using Extents=std::array<size_t,RANK>;


/// @brief Increments a multi-index in row-major (last-index-fastest) order, in-place.
///
/// The iteration pattern is: increment the last index; on overflow, reset it to zero and carry to the previous index.
/// The loop termination guard is `idx[0] != extents[0]` — once the leading index overflows the array is exhausted.
/// The final increment leaves `idx` in an out-of-bounds state, consistent with how unchecked iterators work.
///
/// @note TODO: formulate as a `std::ranges`-compliant range; see `SliceIterator.h` for a draft.
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

/// @brief Rank-0 no-op overload of `incrementMultiIndex` to satisfy template instantiation for degenerate cases.
inline auto& incrementMultiIndex(Extents<0>& idx, Extents<0>) {return idx;}


namespace multiarray {

/// @brief Sentinel type for rank-2 slice notation: `mav(i, _)` returns the i-th row as a rank-1 view.
/// @todo Generalize `operator()` slicing beyond `RANK==2`.
static constexpr struct Placeholder {} _;

} // multiarray


/// @brief A non-owning, strided view over a contiguous buffer. Should be passed by value.
/**
 * Follows pointer-like (reference) semantics regarding constness: a `const MultiArrayView<T,R>`
 * still allows mutation of the underlying data, just as a `T* const` does.
 * All four combinations of `(const) MultiArrayView<(const) T, RANK>` are usable.
 *
 * The flat storage index for multi-index \f$(i_0,\ldots,i_{R-1})\f$ is
 * \f[ \text{offset} + \sum_{\iota=0}^{R-1} s_\iota\, i_\iota \f]
 * where \f$s\f$ are the strides and \f$\text{offset}\f$ accumulates from prior slicing.
 *
 * @todo Alternative design: `MultiArrayBase<StorageType>` parameterized on storage,
 * unifying view and owning cases without inheritance.
 */
template <typename T, size_t RANK>
class MultiArrayView
{
public:
  MultiArrayView(const MultiArrayView&) = default; MultiArrayView(MultiArrayView&&) = default; MultiArrayView() = default;
  MultiArrayView& operator=(const MultiArrayView&) = default; MultiArrayView& operator=(MultiArrayView&&) = default;
  
  /// @brief Constructs a view with explicit extents, strides, offset, and data span. The offset accumulates from prior slicing operations.
  MultiArrayView(Extents<RANK> extents, Extents<RANK> strides, size_t offset, auto&&... dataView)
    : extents{extents}, strides{strides}, offset{offset}, dataView{FWD(dataView)...} {}

  /// @brief Implicit conversion to a const view (`MultiArrayView<const T, RANK>`).
  operator MultiArrayView<const T, RANK>() const requires ( !std::is_const_v<T> ) {return {extents,strides,offset,dataView}; }

  /// @brief Bounds-check all indices against `extents`. Active only in debug builds (`NDEBUG` not defined).
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

  /// @brief Subscript by `Extents<RANK>` multi-index. Always `const` — does not modify the view, only potentially the underlying data.
  /** 
   * The flat index is \f[ o + \sum_{\iota=0}^{R-1} s_\iota\, i_\iota \f]
   * computed via `fold_left` over `zip(idx, strides)`.
   */
  T& operator() (Extents<RANK> idx) const
  {
    std::apply([this] (auto ... args) { checkBounds(args ...); },idx);

    return dataView[std::ranges::fold_left( std::views::zip(idx,strides), offset,
                                            [&] (auto init, auto ids) {return init+get<0>(ids)*get<1>(ids);} ) ];
  }
  
  /// @brief Subscript by individual indices, one per rank. Requires `sizeof...(i) == RANK`.
  T& operator() (std::convertible_to<size_t> auto ... i) const requires (sizeof...(i)==RANK)
  {
    checkBounds(i...);
    size_t idx=0;
    auto s=strides.begin();
    return dataView[ offset + ( ... , [&] { return idx+=(*s++)*i;} () ) ];
  }
  
  /// @brief Rank-1 specialization: single-index subscript.
  T& operator() (std::convertible_to<size_t> auto i) const requires (RANK==1) {checkBounds(i); return dataView[offset+strides[0]*i];}

  /// @brief Rank-2 specialization: two-index subscript.
  T& operator() (std::convertible_to<size_t> auto i, std::convertible_to<size_t> auto j) const requires (RANK==2) {
    checkBounds(i,j); return dataView[offset+strides[0]*i+strides[1]*j];
  }

  /// @brief Rank-2 row slice: `mav(i, _)` returns a rank-1 view of row `i`. @todo Generalize to arbitrary RANK.
  auto operator() (std::convertible_to<size_t> auto i, multiarray::Placeholder) const requires (RANK==2) {
    return MultiArrayView<T,1>{ {extents[1]}, {strides[1]}, offset+i*strides[0], dataView};
  }

  /// @brief Rank-2 column slice: `mav(_, j)` returns a rank-1 view of column `j`. @todo Generalize to arbitrary RANK.
  auto operator() (multiarray::Placeholder, std::convertible_to<size_t> auto j) const requires (RANK==2) {
    return MultiArrayView<T,1>{ {extents[0]}, {strides[0]}, offset+j*strides[1], dataView};
  }

  Extents<RANK> extents, strides;

  size_t offset;
  
  std::span<T> dataView;
  
};


/// @brief `MultiArrayView` is cheap to copy and should always be passed by value.
template <typename T, size_t RANK>
constexpr auto passByValue_v<MultiArrayView<T,RANK>> = true;


/// @brief Convenience alias: a view over immutable elements.
template <typename T, size_t RANK> requires ( !std::is_const_v<T> )
using MultiArrayConstView = MultiArrayView<const T, RANK>;


// MultiArrayView is itself a reference — passing it by value is equivalent to passing by reference.
template <typename T, size_t RANK> requires ( !std::is_const_v<T> )
struct ReferenceMF<MultiArrayView<T,RANK>> : std::type_identity<MultiArrayView<T,RANK>> {};

template <typename T, size_t RANK> requires ( !std::is_const_v<T> )
struct ConstReferenceMF<MultiArrayView<T,RANK>> : std::type_identity<MultiArrayConstView<T,RANK>> {};


/// @brief Debug-mode extent check between two views. Throws `std::runtime_error` on mismatch.
template <typename T1, typename T2, size_t RANK>
void checkExtents(MultiArrayConstView<T1,RANK> m1, MultiArrayConstView<T2,RANK> m2, std::string message)
{
#ifndef   NDEBUG
  if (m1.extents!=m2.extents) throw std::runtime_error("Extent mismatch in "+message+": "+toStringJSON(m1.extents)+" "+toStringJSON(m2.extents));
#endif // NDEBUG
}

/// @brief Element-by-element equality test. Traverses all elements via `incrementMultiIndex`.
/// @note Not spelled `operator==` — this can be an expensive operation.
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

/// @brief Computes row-major (last-index-fastest) strides from extents: `strides[0]=1`, `strides[i]=strides[i-1]*extents[i-1]`.
template <size_t RANK>
auto calculateStrides(Extents<RANK> extents)
{
  Extents<RANK> strides; strides[0]=1;
  for (size_t i=1; i<RANK; ++i) strides[i]=strides[i-1]*extents[i-1];
  return strides;
}

/// @brief Returns the total number of elements: the product of all extents.
template <size_t RANK>
auto calculateExtent(Extents<RANK> extents)
{
  return std::ranges::fold_left(extents,1uz,std::multiplies{});
}

/// @brief Initializer that allocates storage of size `e` without value-initializing elements.
template <typename T>
const auto noInit = [] (size_t e) -> std::vector<T> {return std::vector<T>(e);};

/// @brief Initializer that allocates storage of size `e` and zero-initializes all elements via `T(0.)`.
/// @todo Use `T{}` instead of `T(0.)` for generality beyond numeric types.
template <typename T>
const auto zeroInit = [] (size_t e) -> std::vector<T> {return std::vector<T>(e,T(0.));};


/// @brief Underlying storage type: `std::vector<T>`. May eventually be replaced by `kokkos::View` or `std::valarray`.
template <typename T>
using StorageType = std::vector<T>;


/// @brief Initializer that copies from a sized range, checking that its size matches the requested extent.
template <typename T>
const auto copyInit(const std::ranges::sized_range auto& input) {return [&] (size_t e) {
  if (size(input) != e) throw std::runtime_error("Extent mismatch in MultiArray copyInit: "+std::to_string(size(input))+" "+std::to_string(e));
  return StorageType<T>(begin(input),end(input));
};}


} // multiarray


/// @brief Owning, move-only, rank-`RANK` array. Non-copyable; follows value semantics regarding constness.
/**
 * Derives from `MultiArrayConstView<T,RANK>` to expose the full view interface.
 * The `dataView` span in the base is re-seated to `data_` after each construction or move.
 *
 * @note Known issue: `load()` (Boost.Serialization) deserializes `data_` into a fresh
 * `std::vector` but does not re-seat the inherited `dataView` span — dangling span bug.
 * Move construction is safe since `std::vector` move preserves the buffer address.
 *
 * @note The storage type `std::vector<T>` may eventually be replaced by `kokkos::View`
 * or `std::valarray` for better HPC integration.
 */
template <typename T, size_t RANK>
requires ( !std::is_const_v<T> )
class MultiArray : public MultiArrayConstView<T,RANK>
{
public:
  using StorageType = multiarray::StorageType<T>;
  
  MultiArray(const MultiArray&) = delete; MultiArray& operator=(const MultiArray&) = delete;
  
  MultiArray(MultiArray&&) = default; MultiArray& operator=(MultiArray&&) = default;

  /// @brief Primary constructor. @p initializer is a callable `(size_t totalElements) -> StorageType`.
  MultiArray(Extents<RANK> extents, auto initializer) requires requires (size_t e) { { initializer(e) } -> std::convertible_to<StorageType>;}
    : MultiArrayConstView<T,RANK>{extents,multiarray::calculateStrides(extents),0}, data_{initializer(multiarray::calculateExtent(extents))}
  {
    this->dataView=std::span<T>(data_);
  }

  /// @brief Constructs an uninitialized array with given extents.
  explicit MultiArray(Extents<RANK> extents) : MultiArray{extents,multiarray::noInit<T>} {}

  /// @brief Rank-1 constructor from an existing storage vector (copies the data).
  explicit MultiArray(const StorageType& st) requires (RANK==1) : MultiArray{{st.size()},multiarray::copyInit<T>(st)} {}

  /// @brief Named copy: constructs a new owning array with a deep copy of @p ma's data.
  friend MultiArray copy(const MultiArray& ma) {return MultiArray{ma.extents, [&] (size_t) {return ma.dataStorage();}};}

  /// @brief Element-by-element assignment from a const view of compatible type.
  /// @note Not `operator=` — this is an expensive operation.
  /// @note `assignTo` copies the *source* into `*this` (the name is counterintuitive).
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

  /// @brief Explicit const view, useful where implicit conversion is suppressed.
  MultiArrayConstView<T,RANK> constView() const {return *this;}
  
  /// @brief Explicit mutable view.
  auto mutableView() {return MultiArrayView<T, RANK>{this->extents,this->strides,0,data_};}
  
  /// @brief Implicit conversion to a mutable view.
  operator MultiArrayView<T,RANK>() {return mutableView();}

  /// @brief Non-const element access, forwarding to the base const overload via `const_cast`.
  T& operator()(auto&&... i) {return const_cast<T&>(static_cast<MultiArrayConstView<T,RANK>>(*this)(FWD(i)...)) ;}
  
  /// @brief Const element access. This overload is needed to prevent the non-const overload from shadowing the base.
  const T& operator()(auto&&... i) const {return static_cast<MultiArrayConstView<T,RANK>>(*this)(FWD(i)...) ;}

  /// @brief Direct access to the underlying storage vector.
  StorageType& dataStorage() {return data_;}
  const StorageType& dataStorage() const {return data_;}

  /// @brief Boost.JSON serialization: emits `{"extents": [...], "data": [...]}`.
  friend void tag_invoke( const json::value_from_tag&, json::value& jv, const MultiArray& ma )
  {
    jv = json::object{
      { "extents" , json::value_from(ma.extents) },
      //    { "strides" , ma.strides },
      { "data", json::value_from(ma.data_) }
    };
  }

  /// @brief Boost.JSON deserialization: reconstructs from `{"extents": [...], "data": [...]}`.
  friend MultiArray tag_invoke( const json::value_to_tag< MultiArray >&, const json::value& jv )
  {
    return MultiArray{
      json::value_to<Extents<RANK>>(jv.as_object().at("extents")),
      multiarray::copyInit<T>( json::value_to<StorageType>(jv.as_object().at("data") ) ) };
  }
  
private:
  StorageType data_;

  friend class boost::serialization::access;
  
  /// @brief Boost.Serialization save: archives extents and raw data. @note Transitional — will be superseded by HDF5.
  template<class Archive> void save(Archive& ar, const unsigned int) const {ar & this->extents & data_;}

  /// @brief Boost.Serialization load: restores extents and data, then recomputes strides.
  /// @warning Does not re-seat `dataView` after deserialization — dangling span bug. See class note.
  template<class Archive> void load(Archive& ar, const unsigned int)
  {
    ar & this->extents & data_; this->strides=multiarray::calculateStrides(this->extents); this->offset=0;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
};

/// @brief `MultiArray` owns its data and must not be passed by value.
template <typename T, size_t RANK>
constexpr auto passByValue_v<MultiArray<T,RANK>> = false;


} // cppqedutils