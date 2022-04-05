// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Definition of cppqedutils::SliceIterator together with its helpers}
#ifndef CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED

#include "ArrayTraits.h"
#include "MultiIndexIterator.h"
#include "TMP_Tools.h"

#include <boost/hana.hpp>

#include <boost/range/iterator_range.hpp>

#include <stdexcept>


namespace hana = boost::hana;


namespace cppqedutils {


template <typename RetainedAxes>
constexpr auto hanaSize(RetainedAxes ra=RetainedAxes{}) {return int(hana::size(ra));}


namespace sliceiterator {

using mii::Begin; using mii::End; using tmptools::ordinals;
  
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
void checkConsistency(std::optional<RetainedAxes> ra=std::nullopt) // if the optional argument is not given, RetainedAxes needs to be default-constructible 
{
  using namespace hana; using namespace literals;

  static_assert ( maximum(ra) < int_c<RANK> , "Too large element" );
  static_assert ( minimum(ra) >= 0_c , "Negative element" );

  if constexpr ( Sequence<RetainedAxes>::value ) {
    static constexpr auto sorted = sort(ra);
    
    static_assert ( unique(sorted) == sorted , "Duplicates found" );
  }
}


/// Filters out the indices corresponding to a subsystem
/**
 * \tparamRANK
 * \tparam RetainedAxes compile-time vector specifying the subsystem (cf. \ref specifyingsubsystems)
 * 
 * \param idx the indices to be filtered (of size `RANK`)
 * 
 * \return the indices *not* contained by the subsystem specified by `RetainedAxes`
 * 
 * ** Example **
 *     
 *     #include "SliceIterator.h"
 * 
 *     
 *     int main() {
 *       const IdxTiny<11> idx{21,2,4,10,11,3,5,9,23,22,7};
 *       const auto filteredIdx{filterOut(idx,tmptools::vector_c<3,6,1,9,7>)};
 *       assert( all( filteredIdx==IdxTiny<6>{21,4,11,3,23,7} ) );
 *     }
 * 
 */
template<int RANK, typename RetainedAxes>
auto filterOut(const IdxTiny<RANK>& idx, RetainedAxes ra=RetainedAxes{})
{
  using namespace hana;
  
  int origIndex=0, resIndex=0;
  
  static constexpr auto sorted = [&]() {
    if constexpr ( Sequence<RetainedAxes>::value ) return sort(ra);  
    else return ra;
  }();
  
  IdxTiny<RANK-hanaSize(ra)> res;

  for_each(sorted,[&](auto i) {
    for (; origIndex<i; (++origIndex,++resIndex)) res(resIndex)=idx(origIndex);
    ++origIndex; // skip value found in sorted
  });
  // the final segment:
  for (; origIndex<idx.length(); (++origIndex,++resIndex)) res(resIndex)=idx(origIndex);
  
  return res;

}

/// Helper to transpose
template<int RANK, typename RetainedAxes>
constexpr auto transposingIndeces(RetainedAxes ra=RetainedAxes{})
{
  using namespace hana; using namespace literals;
  
  return fold(
    ordinals<RANK>,make_tuple(tuple_c<int>,0_c),[&](auto state, auto element) {
      if constexpr ( contains(ra,element) )
        return make_tuple(append(state[0_c],ra[state[1_c]]),state[1_c]+1_c);
      else
        return make_tuple(append(state[0_c],element),state[1_c]);
    })[0_c];
}  


/// Performs the “possible permutation” of the retained axes (cf. \ref multiarrayconcept "Synopsis")
/**
 *  This is needed in order that \f$\avr{1,3,6,7,9}\f$ be the corresponding state vector slice, since this is how it is expected in applications. Cf. Composite for further explanations
 * \return Reference to the function argument.
 * 
 * \par Semantics
 * ~~~
 * static const int RANK=11;
 * 
 * ARRAY<RANK> psi;
 * 
 * transpose<RANK>(psi,tmptools::vector_c<3,6,1,9,7>);
 * ~~~
 * is equivalent to ::
 * ~~~
 * psi.transposeSelf(0,3,2,6,4,5,1,9,8,7,10);
 * ~~~
 * that is, in place of the indices specified by the elements of the compile-time vector `RetainedAxes`,
 * the elements of `RetainedAxes` are put, but in the *order* specified by this tuple.
 */
template<int RANK, typename RetainedAxes, typename ARRAY>
ARRAY& transpose(ARRAY& array, RetainedAxes ra=RetainedAxes{})
{
  hana::unpack(transposingIndeces<RANK>(ra),[&](auto && ... i) -> void {array.transposeSelf(std::forward<std::decay_t<decltype(i)>>(i)...);});
  return array;
}


template<int RANK, typename RetainedAxes>
auto slicingTuple(RetainedAxes ra=RetainedAxes{})
{
  using namespace hana; using namespace literals;
  return fold(ordinals<RANK>,hana::make_tuple(hana::make_tuple(),0_c),[&](auto state, auto element) {
    if constexpr ( contains(ra,element) )
      return make_tuple(append(state[0_c],blitz::Range::all()),state[1_c]);
    else
      return make_tuple(append(state[0_c],0),state[1_c]+1_c);
  })[0_c];
}


/// Performs the slicing on an array already transposed
/**
 * \return Reference to resArray
 * 
 * \par Semantics
 * ~~~
 * static const int RANK=11;
 * 
 * ARRAY<RANK> psi;
 * 
 * typedef tmptools::Vector<3,6,1,9,7> Vec;
 * 
 * ResArray<ARRAY,Vec> psiRes;
 * 
 * VecIdxTiny<RANK,Vec> idxTiny;
 * 
 * Indexer<ARRAY,RANK,Vec>::_(psi,psiRes,idxTiny);
 * ~~~
 * is equivalent to
 * ~~~
 * const blitz::Range a(blitz::Range::all());
 * psiRes.reference(psi(0,a,2,a,4,5,a,a,8,a,10));
 * ~~~
 */
template <int RANK, typename RetainedAxes>
class Slicer
{
public:
  static constexpr auto raSize=hanaSize<RetainedAxes>();
  
  template<template <int> class ARRAY>
  static auto& _(const ARRAY<RANK>& array, ///< The array to be sliced
                 ARRAY<raSize>& resArray, ///< The array storing the result.
                 const IdxTiny<RANK-raSize>& idx ///< Set of dummy indices
  )
  {
    using namespace hana; using namespace literals;

    // fill idx_ with the run-time values
    fold(ordinals<RANK>,idx.begin(),[&](auto iter, auto element) {
        if constexpr ( !contains(RetainedAxes{},element) ) idx_[element]=*iter++;
        return iter;
      });
    
    resArray.reference( unpack(idx_,[&](auto && ... i) {return SubscriptMultiArray<ARRAY<RANK>>::_(array,std::forward<std::decay_t<decltype(i)>>(i)...); }) );
    return resArray;      
  }
  
private:
  // preloaded with blitz::Range::all()s at the retained axes
  inline static auto idx_ = slicingTuple<RANK,RetainedAxes>();
                                    
};



/** @} */


////////////////
//
// SliceIterator
//
////////////////

// TODO: Think over: is the following solution based on inheritance to solve the specialization problem optimal? 
// Note: this is NOT allowed (specialization depending on a template parameter)
// template<typename V> SliceIterator<Size_v<V>,V>;

/// Generic worker base for SliceIterator
template<template <int> class ARRAY, int RANK, typename RetainedAxes>
class Base
{
public:
  static constexpr int RANKIDX=RANK-hanaSize<RetainedAxes>();

  template <bool IS_END>
  Base(const ARRAY<RANK>& array, std::bool_constant<IS_END> ie)
    : array_(array), resArray_(), mii_(filterOut<RANK,RetainedAxes>(array.lbound()),filterOut<RANK,RetainedAxes>(array.ubound()),ie)
  {
    if constexpr (!IS_END) {
      array_.reference(array);
      transpose<RANK>(array_,RetainedAxes{});
    }
  }

  void increment() {++mii_;}
  
  auto& operator*() const {return Slicer<RANK,RetainedAxes>::template _<ARRAY>(array_,resArray_,*mii_);}

  friend bool operator==(const Base& i1, const Base& i2) {return i1.mii_==i2.mii_ /* && i1.array_==i2.array_ */;}
  // The user has to ensure that the two arrays are actually the same

  const MultiIndexIterator<RANKIDX>& operator()() const {return mii_;}

private:
  ARRAY<RANK> array_; // By value, so that the necessary transposition specified by RetainedAxes can be readily performed

  mutable ARRAY<hanaSize<RetainedAxes>()> resArray_;
  
  MultiIndexIterator<RANKIDX> mii_;

};


/// used in the case when `size<RetainedAxes> = RANK` and when `RANK = 1`
template<template <int> class ARRAY, typename RetainedAxes>
class BaseTrivial
{
public:
  static constexpr int RANK=hanaSize<RetainedAxes>();

  template <bool IS_END>
  BaseTrivial(const ARRAY<RANK>& array, std::bool_constant<IS_END>)
    : array_(), isEnd_(IS_END)
  {
    if constexpr (!IS_END) array_.reference(array);
  }

  void increment() {if (!isEnd_) isEnd_=true; else throw std::out_of_range("In sliceiterator::BaseTrivial::increment");}

  auto& operator*() const {if (isEnd_) throw std::out_of_range("In sliceiterator::BaseTrivial::operator*"); return array_;}

  friend bool operator==(const BaseTrivial& i1, const BaseTrivial& i2) {return i1.isEnd_==i2.isEnd_;}

protected:
  mutable ARRAY<RANK> array_; // must be mutable in order that the const operator* can return a non-const ARRAY (without mutable here, the auto there resolves to const)

private:
  bool isEnd_;

};


/// used in the case when `size<RetainedAxes> = RANK`
template<template <int> class ARRAY, typename RetainedAxes>
class BaseSpecial : public BaseTrivial<ARRAY,RetainedAxes>
{
public:
  using BaseTrivial<ARRAY,RetainedAxes>::RANK;

  template <bool IS_END>
  BaseSpecial(const ARRAY<RANK>& array, std::bool_constant<IS_END> ie) : BaseTrivial<ARRAY,RetainedAxes>(array,ie)
  {
    if constexpr (!IS_END) transpose<RANK>(this->array_,RetainedAxes{});
  }

};



} // sliceiterator

#define BASE_class std::conditional_t<RANK==1,\
                                      sliceiterator::BaseTrivial<ARRAY,RetainedAxes>,\
                                      std::conditional_t<RANK==hanaSize<RetainedAxes>(),\
                                                         sliceiterator::BaseSpecial<ARRAY,RetainedAxes>,\
                                                         sliceiterator::Base<ARRAY,RANK,RetainedAxes>>>


/// SliceIterator
/**
 * \tparam ARRAY the array to operate on, practically either a blitz::Array, or a quantumdata::StateVector or quantumdata::DensityOperator
 * \tparam RANK positive integer standing for the number of elementary Hilbert spaces
 * \tparam RetainedAxes compile-time tuple holding the *retained index positions* like \f$\avr{3,6,1,9,7}\f$ \ref retainedindexpositionsdefined "here". (Cf. \ref specifyingsubsystems)
 * 
 * To understand the template parameters, cf. also \ref multiarrayconcept.
 * 
 * Model of [ForwardIterator](http://www.cplusplus.com/reference/std/iterator/ForwardIterator/).
 * 
 * This iterator is implemented in terms of a MultiIndexIterator, and hence it can be initialized to either the beginning or the end of the sequence of retained indices
 * 
 * This class is at the absolute heart of the framework as it is indispensable to implement 
 * - composite quantum systems (Composite and BinarySystem)
 * - iterations over (multi)matrix rows and columns to get (multi)vectors
 * - diagonal iteration of LazyDensityOperator
 * 
 * This said, it is never really used directly in the framework, but rather through the maker functions
 * sliceiterator::begin, sliceiterator::end, and sliceiterator::fullRange in standard or
 * \refBoost{Boost.Range algorithms,range/doc/html/range/reference/algorithms.html}.
 * 
 * Quite generally, by iterating through all the combinations of indices *not* belonging to the given subsystem (dummy indices) and when dereferenced returning the corresponding slice,
 * it can be used to implement the action of operators in extended (and/or permutated) Hilbert spaces.
 * 
 * \Semantics
 * 
 * Sticking to the example \ref retainedindexpositionsdefined "here", assume that the function
 * ~~~
 * void actWithA(CArray<5>&);
 * ~~~
 * implements the action of the operator \f$A\f$ on a state vector of rank 5. Then the action on the extended Hilbert space can be calculated as
 * ~~~
 * void actOnExtended(CArray<11>& psi)
 * {
 *   for (auto&& psiS : cppqedutils::sliceiterator::fullRange<tmptools::Vector<3,6,1,9,7> >(psi)) actWithA(psiS);
 * }
 * ~~~
 * 
 * \see cppqedutils::sliceiterator::fullRange
 * 
 * \todo Implement a default version of SliceIterator for the case when neither slicing nor transposition is necessary, that is when `V` is equivalent to a range<0,RANK-1>.
 * This will require further compile-time implementation selection.
 *
 */
template<template <int> class ARRAY, int RANK, typename RetainedAxes>
class SliceIterator 
  : public boost::forward_iterator_helper<SliceIterator<ARRAY,RANK,RetainedAxes>,ARRAY<hanaSize<RetainedAxes>()>>, // The inheritance has to be public here because of needed inherited types 
    public BASE_class
{
public:
  typedef BASE_class Base;

#undef  BASE_class

  SliceIterator& operator++() {Base::increment(); return *this;} ///< For the ForwardIterator concept

  /// Can be initialized either to the beginning or the end of the sequence of dummy-index combinations
  /** \tparam IS_END governs the end-ness */
  template<bool IS_END>
  SliceIterator(const ARRAY<RANK>& array, std::bool_constant<IS_END> ie) : Base(array,ie) {}
  
private:
  static void checkConsistency() {sliceiterator::checkConsistency<RANK,RetainedAxes>();}

};


namespace sliceiterator {


/// Iterator to the beginning of the sequence \related SliceIterator – template parameters have the same meaning as there \return SliceIterator
template<typename RetainedAxes, template <int> class ARRAY, int RANK>
auto begin(const ARRAY<RANK>& array) {return SliceIterator<ARRAY,RANK,RetainedAxes>(array,Begin());}

/// Iterator to the end of the sequence \related SliceIterator – template parameters have the same meaning as there \return SliceIterator
template<typename RetainedAxes, template <int> class ARRAY, int RANK>
auto end(const ARRAY<RANK>& array) {return SliceIterator<ARRAY,RANK,RetainedAxes>(array,End());}

/// \refBoost{Boost.Range,range/doc/html/range/reference/utilities/iterator_range.html}-compliant full range of slice iterators \related SliceIterator
/**
 * It corresponds to all the possible combinations of dummy indices (\f$\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},\f$…) \ref retainedindexpositionsdefined "here".
 * \tparam RetainedAxes the retained index positions \tparam A array type taken itself as a parameter to ease template-parameter inference 
 */
template<typename RetainedAxes, template <int> class ARRAY, int RANK>
auto fullRange(const ARRAY<RANK>& array)
{
  return boost::iterator_range<SliceIterator<ARRAY,RANK,RetainedAxes>>(begin<RetainedAxes,ARRAY>(array),end<RetainedAxes,ARRAY>(array));
/*                               SliceIterator<ARRAY,RANK,V>,
                                 std::ranges::subrange_kind::unsized>(begin<V,ARRAY>(array),end<V,ARRAY>(array));*/
}


/** \page specifyingsubsystems Specifying subsystems
 * 
 * It is assumed throughout that RetainedAxes is given as a `hana::tuple_c<int,...>` or `hana::range_c<int,...>`
 * Since the type of such an object is a `final` class, it is safe to assume that it is default-constructible 
 * 
 * Many constructs in the framework require the specification of a subsystem of a quantum system of arbitrary arity @ compile time.
 * The main example is retained index positions for slicing (cf. \ref multiarrayconcept).
 * It is the template parameter `RetainedAxes`, a compile-time sequence, that specifies the subsystem.
 * 
 * \par Preconditions
 * 
 * `RetainedAxes` must not contain
 * - negative values,
 * - values not smaller than the arity, and
 * - duplicate values.
 * These imply that the size of `RetainedAxes` cannot be greater than `RANK`.
 * 
 * These are checked for @ compile time, and any violation is signalled by assertions from Boost.Hana
 * 
 * \see checkConsistency
 * 
 */


} // sliceiterator

} // cppqedutils

#endif // CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED
