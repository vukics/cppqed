// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Definition of cppqedutils::SliceIterator together with its helpers}
#ifndef CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_SLICEITERATOR_H_INCLUDED

#include "MultiIndexIterator.h"
#include "TMP_Tools.h"

#include <boost/range.hpp>

#include <boost/mpl/size.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/max_element.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>

#include <boost/fusion/mpl.hpp> // include all, otherwise difficult to find out what is actually needed
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/algorithm/iteration/fold.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>

#include <stdexcept>


namespace cppqedutils {


namespace sliceiterator {
  
/// A forwarding metafunction to boost::mpl::size
template <typename V>
constexpr size_t Size_v = boost::mpl::size<V>::value;

template <template <int> class ARRAY, typename V> using ResArray=ARRAY<Size_v<V>>;

template <int RANK, typename V> using VecIdxTiny=IdxTiny<RANK-Size_v<V>>; 
// note that Array::lbound and ubound return a TinyVector<int,...>, so that here we have to use int as well.

using mii::Begin; using mii::End;


/// Checking the consistency of template arguments for use in slicing
/**
 * - size of `V` must not be larger than `RANK`
 * - `V` must not contain duplicated elements
 * - all elements of `V` must be smaller than `RANK`
 * 
 * \see \ref specifyingsubsystems
 * 
 */
template<int RANK, typename V> struct ConsistencyChecker
{
  static_assert( RANK >= Size_v<V> , "Indexer with nonpositive RANK." );

  typedef typename boost::mpl::sort<V>::type SortedV;
  static_assert( boost::mpl::equal<typename boost::mpl::unique<SortedV,boost::is_same<boost::mpl::_1,boost::mpl::_2> >::type,SortedV>::value ,
                 "cppqedutils::SliceIterator inconsistent vector" );
  static_assert( boost::mpl::deref<typename boost::mpl::max_element<V>::type>::type::value < RANK , 
                 "cppqedutils::SliceIterator vector out of range" );
};


/// Specialization necessary since boost::mpl::range_c cannot be sorted
template<int RANK, int I1, int I2>  struct ConsistencyChecker<RANK,boost::mpl::range_c<int,I1,I2> >
{
  static_assert( I1>=0 && I2<RANK , 
                 "cppqedutils::SliceIterator vector out of range" );
};

/** \cond SPECIALIZATION */
template<int RANK, int N, int Nbeg> struct ConsistencyChecker<RANK,tmptools::Range<N,Nbeg> > : private ConsistencyChecker<RANK,boost::mpl::range_c<int,Nbeg,Nbeg+N> > {};
template<int RANK, int N>           struct ConsistencyChecker<RANK,tmptools::Ordinals<N>   > : private ConsistencyChecker<RANK,tmptools::  Range  <N,0>             > {};
/** \endcond */


/// Filters out the indices corresponding to a subsystem
/**
 * \tparamRANK
 * \tparam V compile-time vector specifying the subsystem (cf. \ref specifyingsubsystems)
 * 
 * \param idx the indices to be filtered (of number `RANK`)
 * 
 * \return the indices *not* contained by the subsystem specified by `V` (of number `RANK-Size_v<V>`)
 * 
 */
template<int RANK, typename V>
auto filterOut(const IdxTiny<RANK>& v)
{
  namespace mpl=boost::mpl;

  unsigned curr=0;
  VecIdxTiny<RANK,V> res;

  auto helper=[&](auto i){res(curr++)=v(decltype(i)::value);};
  
  mpl::for_each<mpl::filter_view<tmptools::Ordinals<RANK>,mpl::not_<tmptools::numerical_contains<V,mpl::_> > > >(helper);
  
  return res;

}


/**
 * \defgroup helperstoiterator Helpers to BlitzArraySliceIterator
 * 
 * Sometimes used on their own as well, so they are exposed in the header file.
 * 
 * \tparamRANK \tparamV
 * 
 * \pre `Size_v<V> <= RANK`
 * 
 * The technique of using (non-templated) static worker functions of class templates is meant to allow partial specializations, which are not possible for function templates.
 * 
 * \internal These are just rudimentary definitions, the actual definitions being partial specializations of Transposer::_ and Indexer::_ along RANK 
 * (cf. trailing part of `BlitzArraySliceIterator.tcc`) The functions will throw if a partial specialization for the given RANK does not exist. \endinternal
 * 
 * \see \ref iteratorimplementation,*this
 * 
 * @{
 */

/// Class performing the “possible permutation” of the retained indices (cf. \ref multiarrayconcept "Synopsis").
template<template <int> class ARRAY, int RANK, typename V>
class Transposer
{
public:
  /// Static worker
  /**
   * Transposition corresponding to the "possible permutation" of the retained indices (cf. \ref multiarrayconcept "Synopsis"), which is necessary in order that \f$\avr{1,3,6,7,9}\f$
   * be the corresponding state vector slice, since this is how it is expected in applications. Cf. Composite for further explanations.
   * 
   * \return Simply the reference to the function argument.,*this
   * 
   * \par Semantics
   * ~~~
   * static const int RANK=11;
   * 
   * ARRAY<RANK> psi;
   * 
   * typedef tmptools::Vector<3,6,1,9,7> Vec;
   * 
   * Transposer<ARRAY,RANK,Vec>::_(psi);
   * ~~~
   * is equivalent to ::
   * ~~~
   * psi.transposeSelf(0,3,2,6,4,5,1,9,8,7,10);
   * ~~~
   * that is, in place of the indices specified by the elements of the compile-time vector `Vec`, the elements of `Vec` are put, but in the *order* specified by `Vec`.
   */
  static ARRAY<RANK>& _(ARRAY<RANK>&);

};

/// Performs the slicing on an array already transposed by Transposer.
template<template <int> class ARRAY, int RANK, typename V>
class Indexer
{
public:
  /// Static worker
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
  static ResArray<ARRAY,V>& _(const ARRAY<RANK>&, ///< The array to be sliced
                              ResArray<ARRAY,V>&, ///< The array storing the result. Necessary for the same reason why SliceIterator::operator* returns a reference
                              const VecIdxTiny<RANK,V>& ///< Set of dummy indices
                             );

};


/** @} */


//////////////////////////
//
// BlitzArraySliceIterator
//
//////////////////////////

// TODO: Think over: is the following solution based on inheritance to solve the specialization problem optimal?  Note: this is NOT allowed (specialization depending on a template parameter)
// template<typename V> SliceIterator<Size_v<V>,V>;

template<template <int> class, typename>
class BaseTrivial;

template<template <int> class, typename>
class BaseSpecial;

template<template <int> class, int, typename>
class Base;

} // sliceiterator

#define BASE_class std::conditional_t<RANK==1,\
                                      sliceiterator::BaseTrivial<ARRAY,V>,\
                                      std::conditional_t<RANK==sliceiterator::Size_v<V>,\
                                                         sliceiterator::BaseSpecial<ARRAY,V>,\
                                                         sliceiterator::Base<ARRAY,RANK,V>>>


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
template<template <int> class ARRAY, int RANK, typename V>
class SliceIterator 
  : public boost::forward_iterator_helper<SliceIterator<ARRAY,RANK,V>,sliceiterator::ResArray<ARRAY,V>>, // The inheritance has to be public here because of needed inherited types 
    public BASE_class,
    private sliceiterator::ConsistencyChecker<RANK,V>
{
public:
  typedef BASE_class Base;

#undef  BASE_class

  SliceIterator& operator++() {Base::increment(); return *this;} ///< For the ForwardIterator concept

  /// Can be initialized either to the beginning or the end of the sequence of dummy-index combinations
  /** \tparam IS_END governs the end-ness */
  template<bool IS_END>
  SliceIterator(const ARRAY<RANK>& array, std::bool_constant<IS_END> ie) : Base(array,ie) {}

};


namespace sliceiterator {


/// Iterator to the beginning of the sequence \related SliceIterator – template parameters have the same meaning as there \return SliceIterator
template<typename V, template <int> class ARRAY, int RANK>
auto begin(const ARRAY<RANK>& array) {return SliceIterator<ARRAY,RANK,V>(array,Begin());}

/// Iterator to the end of the sequence \related SliceIterator – template parameters have the same meaning as there \return SliceIterator
template<typename V, template <int> class ARRAY, int RANK>
auto end(const ARRAY<RANK>& array) {return SliceIterator<ARRAY,RANK,V>(array,End());}

/// \refBoost{Boost.Range,range/doc/html/range/reference/utilities/iterator_range.html}-compliant full range of slice iterators \related SliceIterator
/**
 * It corresponds to all the possible combinations of dummy indices (\f$\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},\f$…) \ref retainedindexpositionsdefined "here".
 * \tparam V the retained index positions \tparam A array type taken itself as a parameter to ease template-parameter inference 
 */
template<typename V, template <int> class ARRAY, int RANK>
auto fullRange(const ARRAY<RANK>& array)
{
  return boost::iterator_range<SliceIterator<ARRAY,RANK,V>>(begin<V,ARRAY>(array),end<V,ARRAY>(array));
}


/// Generic worker base for SliceIterator
template<template <int> class ARRAY, int RANK, typename V>
class Base
{
public:
  static const int RANKIDX=RANK-Size_v<V>;

  template <bool IS_END>
  Base(const ARRAY<RANK>& array, std::bool_constant<IS_END> ie)
    : array_(array), resArray_(), mii_(filterOut<RANK,V>(array.lbound()),filterOut<RANK,V>(array.ubound()),ie)
  {
    if constexpr (!IS_END) {
      array_.reference(array);
      Transposer<ARRAY,RANK,V>::_(array_);
    }
  }

  void increment() {++mii_;}
  
  auto& operator*() const {return Indexer<ARRAY,RANK,V>::_(array_,resArray_,*mii_);}

  friend bool operator==(const Base& i1, const Base& i2) {return i1.mii_==i2.mii_ /* && i1.array_==i2.array_ */;}
  // The user has to ensure that the two arrays are actually the same

  const MultiIndexIterator<RANKIDX>& operator()() const {return mii_;}

private:
  ARRAY<RANK> array_; // By value, so that the necessary transposition specified by V can be readily performed

  mutable ResArray<ARRAY,V> resArray_;
  
  MultiIndexIterator<RANKIDX> mii_;

};


/// used in the case when `Size_v<V> = RANK` and when `RANK = 1`
template<template <int> class ARRAY, typename V>
class BaseTrivial
{
public:
  static const int RANK=Size_v<V>;

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


/// used in the case when `Size_v<V> = RANK`
template<template <int> class ARRAY, typename V>
class BaseSpecial : public BaseTrivial<ARRAY,V>
{
public:
  using BaseTrivial<ARRAY,V>::RANK;

  template <bool IS_END>
  BaseSpecial(const ARRAY<RANK>& array, std::bool_constant<IS_END> ie) : BaseTrivial<ARRAY,V>(array,ie)
  {
    if constexpr (!IS_END) Transposer<ARRAY,Size_v<V>,V>::_(this->array_);
  }

};


/** \page iteratorimplementation Notes on the implementation of SliceIterator
 * 
 * \tableofcontents 
 * 
 * Transposer and Indexer are implemented in such a way that a partial template specialization is provided for each possible `RANK` (up to #BLITZ_ARRAY_LARGEST_RANK), 
 * and the corresponding code is automatically generated by the \refBoost{Boost.Preprocessor,preprocessor/doc/index.html} library.
 * This can be seen in the trailing part of BlitzArraySliceIterator.tcc. To actually see what code is generated, this file needs to be preprocessed.
 * Issue the following command from the root directory of the distribution:
 * ~~~{.sh}
 * g++ -P -E -Iutils/ utils/BlitzArraySliceIterator.tcc | tail -n286
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
 * This is done by the following snippet in `utils/BlitzArraySliceIterator.tcc`: \dontinclude BlitzArraySliceIterator.tcc
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
