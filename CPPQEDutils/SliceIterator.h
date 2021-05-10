// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Definition of cppqedutils::SliceIterator together with its helpers}
#ifndef CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED

#include "ArrayTraits.h"
#include "MultiIndexIterator.h"
#include "TMP_Tools.h"

#include <boost/range.hpp>

#include <boost/hana.hpp>

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
void checkConsistency(RetainedAxes ra=RetainedAxes{})
{
  using namespace hana; using namespace literals;

  BOOST_HANA_CONSTANT_ASSERT_MSG ( maximum(ra) < int_c<RANK> , "Too large element" );
  BOOST_HANA_CONSTANT_ASSERT_MSG ( minimum(ra) >= 0_c , "Negative element" );

  if constexpr ( Sequence<RetainedAxes>::value ) {
    static constexpr auto sorted = sort(ra);
    
    BOOST_HANA_CONSTANT_ASSERT_MSG ( unique(sorted) == sorted , "Duplicates found" );
  }
}


/// Filters out the indices corresponding to a subsystem
/**
 * \tparamRANK
 * \tparam RetainedAxes compile-time vector specifying the subsystem (cf. \ref specifyingsubsystems)
 * 
 * \param idx the indices to be filtered (of number `RANK`)
 * 
 * \return the indices *not* contained by the subsystem specified by `RetainedAxes`
 * 
 * ** Example **
 *     
 *     #include "SliceIterator.h"
 * 
 *     const int RANK=11;
 *     constexpr auto v_c=vector_c<3,6,1,9,7>;
 *     
 *     int main() {
 *       const IdxTiny<RANK> idx{21,2,4,10,11,3,5,9,23,22,7};
 *       const auto filteredIdx{filterOut(idx,v_c)};
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
 * typedef tmptools::Vector<3,6,1,9,7> Vec;
 * 
 * Transposer<ARRAY<RANK>,Vec>::_(psi);
 * ~~~
 * is equivalent to ::
 * ~~~
 * psi.transposeSelf(0,3,2,6,4,5,1,9,8,7,10);
 * ~~~
 * that is, in place of the indices specified by the elements of the compile-time vector `Vec`, the elements of `Vec` are put, but in the *order* specified by `Vec`.
 */
template<typename RetainedAxes, typename ARRAY>
ARRAY& transpose(ARRAY& array, RetainedAxes ra=RetainedAxes{})
{
  hana::unpack(transposingIndeces<Rank_v<ARRAY>>(ra),[&](auto && ... i) -> void {array.transposeSelf(std::forward<std::decay_t<decltype(i)>>(i)...);});
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
 * \return Reference to resArray,*this
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
      transpose(array_,RetainedAxes{});
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
    if constexpr (!IS_END) transpose(this->array_,RetainedAxes{});
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
 * \tparam V compile-time vector holding the *retained index positions* like \f$\avr{3,6,1,9,7}\f$ \ref retainedindexpositionsdefined "here". (Cf. \ref specifyingsubsystems)
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
 *   boost::for_each(cppqedutils::sliceiterator::fullRange<tmptools::Vector<3,6,1,9,7> >(psi),
 *                   actWithA);
 * }
 * ~~~
 * 
 * \see cppqedutils::sliceiterator::fullRange and the \refBoostConstruct{for_each,range/doc/html/range/reference/algorithms/non_mutating/for_each.html} algorithm of Boost.Range.
 * 
 * For further basic examples of usage cf. `utils/testsuite/BlitzArraySliceIterator.cc` & `utils/testsuite/BlitzArraySliceIteratorTMP.cc`.
 * 
 * \see \ref iteratorimplementation
 * 
 * \todo Implement a default version of SliceIterator for the case when neither slicing nor transposition is necessary, that is when `V` is equivalent to a range<0,RANK-1>.
 * This will require further compile-time implementation selection.
 *
 * \todo Refine the iterator category according to the \refBoost{New-Style iterator concepts,iterator/doc/index.html#new-style-iterators}.
 * The point is that a multi-array is not a container of slices, so SliceIterator is definitely not a standard iterator. It seems rather like a proxy iterator.
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
}



/** \page iteratorimplementation Notes on the implementation of SliceIterator
 * 
 * \tableofcontents 
 * 
 * Transposer and Indexer are implemented in such a way that a partial template specialization is provided for each possible `RANK` (up to #BLITZ_ARRAY_LARGEST_RANK), 
 * and the corresponding code is automatically generated by the \refBoost{Boost.Preprocessor,preprocessor/doc/index.html} library.
 * This can be seen in the trailing part of BlitzArraySliceIterator.h. To actually see what code is generated, this file needs to be preprocessed.
 * Issue the following command from the root directory of the distribution:
 * ~~~{.sh}
 * g++ -P -E -Iutils/ utils/BlitzArraySliceIterator.h | tail -n286
 * ~~~
 * To store and manipulate the heterogenous collection of slicing indices like `0,a,2,a,4,5,a,a,8,a,10` in Indexer::_, the `vector` class of the 
 * \refBoost{Boost.Fusion,fusion/doc/html/index.html} library is used.
 * 
 * SliceIterator is implemented in terms of the above two helper classes. Each SliceIterator invokes a Transposer::_ at its construction, and – if `RANK` is larger
 * than `Size_v<V>` – an Indexer::_ @ the point of its dereferencing when the actual slicing occurs.
 * 
 * SliceIterator is a forward iterator, implemented with the help of \refBoostConstruct{forward_iterator_helper,utility/operators.htm#iterator} from Boost.Operator.
 * For this to work, we need to define only 3 operations:
 * -# Comparison for equality
 * -# Prefix increment
 * -# Dereferencing
 * A special implementation is needed when the size of the compile-time vector `V` equals `RANK` because in this case actually no slicing takes place,
 * only transposition. For this, as at several other places in the framework, we apply conditional inheritance: SliceIterator inherits from either of two classes 
 * (details::Base or details::BaseSpecial).
 * 
 * The iteration over dummy indices is implemented with the help of cppqedutils::MultiIndexIterator.
 * 
 * A metaprogramming example {#metaprogrammingexample}
 * =========================
 * 
 * In the following we analyse a metaprogramming example typical for the framework: how the compile-time vector `0,3,2,6,4,5,1,9,8,7,10`
 * for the self-transposition in Transposer::_ is prepared.
 * 
 * This is done by the following snippet in `utils/BlitzArraySliceIterator.h`: \dontinclude BlitzArraySliceIterator.h
 * \skip namespace namehider {
 * \until fold
 * We are using the \refBoostConstruct{fold,mpl/doc/refmanual/fold.html} metaalgorithm from Boost.MPL.
 * Here, it iterates over the sequence of ordinals (tmptools::Ordinals) between `0` and `RANK-1`.
 * \until vector_c
 * The initial state for the fold algorithm is an empty \refBoost{compile-time vector of integers,mpl/doc/refmanual/vector-c.html}
 * and the \refBoost{iterator,mpl/doc/refmanual/begin.html} pointing to the first element of the compile-time vector `V`. 
 * These two are “zipped” into a \refBoost{compile-time pair,mpl/doc/refmanual/pair.html}.
 * At the end, the first element of this pair will hold the result.
 * 
 * The rest expresses the forward operation of the fold algorithm in the form of a \refBoost{compile-time lambda expression,mpl/doc/tutorial/handling-placeholders.html}.
 * After each step, the new state will again be a pair composed of
 * \until _2
 * \until _1
 * \until _2
 * \until >,
 * i. the vector is \refBoost{augmented,mpl/doc/refmanual/push-back.html} either by the current element in `V` pointed to by the iterator (Line 13), or the current ordinal (Line 14),
 * depending on whether `V` tmptools::numerical_contains this ordinal.
 * \until second
 * \until >
 * \until >
 * ii. the iterator is \refBoost{advanced,mpl/doc/refmanual/next.html} (Line 18) if the same numerical containment criterion is met,
 * otherwise it is left untouched for the same element to be considered again (Line 19).
 * \note In template metaprogramming, “calling” tmptools::numerical_contains<V,mpl::_2> twice is no waste because the second time no further template 
 * instantiations are needed, the compiler will use the ones already instantiated the first time.
 * 
 * \until TransposerMeta
 * \until {};
 * the \refBoost{first element,mpl/doc/refmanual/trivial-metafunctions-summary.html#first} of the resulting pair of the above algorithm is picked. As it is easy to verify,
 * ~~~
 * TransposerMeta<11,tmptools::Vector<3,6,1,9,7> >::type
 * ~~~
 * will be an
 * ~~~
 * mpl::vector_c<int,0,3,2,6,4,5,1,9,8,7,10>
 * ~~~  
 */

/** \page specifyingsubsystems Specifying subsystems
 * 
 * Many constructs in the framework require the specification of a subsystem of a multiary quantum system @ compile time. The main example is retained index positions for slicing
 * (cf. \ref multiarrayconcept). It is the template parameter `V`, a compile-time sequence, that specifies the subsystem.
 * 
 * \par Example models
 * 
 * tmptools::Vector and \refBoostConstruct{range_c,mpl/doc/refmanual/range-c.html} from Boost.MPL.
 * 
 * \par Preconditions
 * 
 * - `Size_v<V>` must not be larger than the arity
 * - `V` must not “contain”
 *   - negative values,
 *   - values not smaller than the arity, and
 *   - duplicate values.
 * 
 * These are checked for @ compile time, and any violation is signalled by more or less intelligent compiler errors generated with the help of
 * \refBoost{Boost.MPL's static assertions,mpl/doc/refmanual/asserts.html}.
 * 
 * \see ConsistencyChecker, tmptools::Vector
 * 
 */


} // sliceiterator

} // cppqedutils

#endif // CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED
