// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFile{Definition of blitzplusplus::basi::Iterator together with its helpers}
#ifndef CPPQEDCORE_UTILS_BLITZARRAYSLICEITERATOR_H_INCLUDED
#define CPPQEDCORE_UTILS_BLITZARRAYSLICEITERATOR_H_INCLUDED

#include "BlitzArraySliceIteratorFwd.h"

#include "BlitzArrayTraits.h"
#include "MultiIndexIterator.h"
#include "TMP_Tools.h"

#include "Exception.h"

#include <boost/mpl/size.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/max_element.hpp>

#include <boost/range.hpp>

#include <list>



namespace blitzplusplus {


using cpputils::Rank;


/// The name of the namespace stands for <strong>B</strong>litz<strong>A</strong>rray<strong>S</strong>lice<strong>I</strong>terator
namespace basi {

/// A forwarding metafunction to boost::mpl::size
template <typename V> struct Size : boost::mpl::size<V> {};

/// Contains template aliases for use in BlitzArraySliceIterator.h & BlitzArraySliceIterator.tcc
namespace ttd {

template <int RANK, bool IS_CONST> using ConditionalConstCArray=typename tmptools::ConditionalAddConst<CArray<RANK>,IS_CONST>::type;

template <typename V> using ResCArray=CArray<Size<V>::value>;

template <typename V, bool IS_CONST> using ConditionalConstResCArray=ConditionalConstCArray<Size<V>::value,IS_CONST>;

template <typename I, typename V, bool IS_CONST> using ForwardIteratorHelper=boost::forward_iterator_helper<I,ConditionalConstResCArray<V,IS_CONST> >;

template <int RANK, typename V> using VecIdxTiny=IdxTiny<RANK-Size<V>::value>; 
// note that Array::lbound and ubound return a TinyVector<int,...>, so that here we have to use int as well.

}

using cpputils::mii::Begin; using cpputils::mii::End;


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
  static_assert( RANK >= Size<V>::value , "Indexer with nonpositive RANK." );

  typedef typename boost::mpl::sort<V>::type SortedV;
  static_assert( boost::mpl::equal<typename boost::mpl::unique<SortedV,boost::is_same<boost::mpl::_1,boost::mpl::_2> >::type,SortedV>::value , "basi::Iterator inconsistent vector" );
  static_assert( boost::mpl::deref<typename boost::mpl::max_element<V>::type>::type::value < RANK , "basi::Iterator vector out of range" );
};


/// Specialization necessary since boost::mpl::range_c cannot be sorted
template<int RANK, int I1, int I2>  struct ConsistencyChecker<RANK,boost::mpl::range_c<int,I1,I2> > { static_assert( I1>=0 && I2<RANK , "basi::Iterator vector out of range" ); };

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
 * \return the indices *not* contained by the subsystem specified by `V` (of number `RANK-Size<V>::value`)
 * 
 */
template<int RANK, typename V>
const ttd::VecIdxTiny<RANK,V>
filterOut(const IdxTiny<RANK>& idx);



/////////////////////////////////
//
// The essential part of the work
//
/////////////////////////////////


/**
 * \defgroup helperstoiterator Helpers to BlitzArraySliceIterator
 * 
 * Sometimes used on their own as well, so they are exposed in the header file.
 * 
 * \tparamRANK \tparamV
 * 
 * \pre `Size<V> <= RANK`
 * 
 * The technique of using (non-templated) static worker functions of class templates is meant to allow partial specializations, which are not possible for function templates.
 * 
 * \internal These are just rudimentary definitions, the actual definitions being partial specializations of Transposer::transpose and Indexer::index along RANK 
 * (cf. trailing part of `BlitzArraySliceIterator.tcc`) The functions will throw if a partial specialization for the given RANK does not exist. \endinternal
 * 
 * \see \ref iteratorimplementation
 * 
 * @{
 */

/// Exception thrown if the partial specialization of Transposer or Indexer for the given RANK does not exist.
template<int RANK>
class TransposerOrIndexerRankTooHighException : public cpputils::Exception {};

/// Class performing the “possible permutation” of the retained indices (cf. \ref multiarrayconcept "Synopsis").
template<int RANK, typename V>
class Transposer
{
public:
  /// Static worker
  /**
   * Transposition corresponding to the "possible permutation" of the retained indices (cf. \ref multiarrayconcept "Synopsis"), which is necessary in order that \f$\avr{1,3,6,7,9}\f$
   * be the corresponding state vector slice, since this is how it is expected in applications. Cf. Composite for further explanations.
   * 
   * \return Simply the reference to the function argument.
   * 
   * \par Semantics
   * ~~~
   * static const int RANK=11;
   * 
   * CArray<RANK> psi;
   * 
   * typedef tmptools::Vector<3,6,1,9,7> Vec;
   * 
   * Transposer<RANK,Vec>::transpose(psi);
   * ~~~
   * is equivalent to ::
   * ~~~
   * psi.transposeSelf(0,3,2,6,4,5,1,9,8,7,10);
   * ~~~
   * that is, in place of the indices specified by the elements of the compile-time vector `Vec`, the elements of `Vec` are put, but in the *order* specified by `Vec`.
   */
  static 
  CArray<RANK>&
  transpose(CArray<RANK>&)
  {throw TransposerOrIndexerRankTooHighException<RANK>();}

};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

/// Performs the slicing on an array already transposed by Transposer.
template<int RANK, typename V>
class Indexer : public Transposer<RANK,V>
{
public:
  /// Static worker
  /**
   * \return Reference to resArray
   * 
   * \par Semantics
   * ~~~
   * static const int RANK=11;
   * 
   * CArray<RANK> psi;
   * 
   * typedef tmptools::Vector<3,6,1,9,7> Vec;
   * 
   * ttd::ResCArray<Vec> psiRes;
   * 
   * ttd::VecIdxTiny<RANK,Vec> idxTiny;
   * 
   * Indexer<RANK,Vec>::index(psi,psiRes,idxTiny);
   * ~~~
   * is equivalent to
   * ~~~
   * const blitz::Range a(blitz::Range::all());
   * psiRes.reference(psi(0,a,2,a,4,5,a,a,8,a,10));
   * ~~~
   */
  static
  ttd::ResCArray<V>&
  index(CArray<RANK>& array, ///< The array to be sliced
        ttd::ResCArray<V>& resArray, ///< The array storing the result
        const ttd::VecIdxTiny<RANK,V>& idx ///< Set of dummy indices
       )
  {throw TransposerOrIndexerRankTooHighException<RANK>();}

};

#pragma GCC diagnostic pop

/** @} */


//////////////////////////
//
// BlitzArraySliceIterator
//
//////////////////////////

// Think over: is the following solution based on inheritance to solve the specialization problem optimal?  Note: this is NOT allowed (specialization depending on a template parameter)
// template<typename V, bool IS_CONST> Iterator<Size<V>::value,V,IS_CONST>;

namespace details {

template<typename, bool>
class BaseTrivial;

template<typename, bool>
class BaseSpecial;

template<int, typename, bool>
class Base;

} // details

#define BASE_class boost::mpl::if_c<RANK==1,\
                                    details::BaseTrivial<V,IS_CONST>,\
                                    typename boost::mpl::if_c<RANK==Size<V>::value,\
                                                              details::BaseSpecial<V,IS_CONST>,\
                                                              details::Base<RANK,V,IS_CONST>\
                                                             >::type\
                                   >::type

/// BlitzArraySliceIterator
/**
 * \tparam RANK positive integer standing for the number of elementary Hilbert spaces
 * \tparam V compile-time vector holding the *retained index positions* like \f$\avr{3,6,1,9,7}\f$ \ref retainedindexpositionsdefined "here". (Cf. \ref specifyingsubsystems)
 * \tparam IS_CONST governs the constness of the class
 * 
 * To understand the template parameters, cf. also \ref multiarrayconcept.
 * 
 * Model of [ForwardIterator](http://www.cplusplus.com/reference/std/iterator/ForwardIterator/). Can be both const and non-const iterator depending on the last template argument.
 * 
 * This iterator is implemented in terms of a cpputils::MultiIndexIterator, and hence it can be initialized to either the beginning or the end of the “sequence”.
 * 
 * This class is at the absolute heart of the framework as it is indispensable to implement 
 * - composite quantum systems (Composite and BinarySystem)
 * - iterations over (multi)matrix rows and columns to get (multi)vectors
 * - quantumdata::ldo::DiagonalIterator
 * 
 * This said, it is never really used directly in the framework, but rather through the maker functions
 * blitzplusplus::basi::begin, blitzplusplus::basi::end, and blitzplusplus::basi::fullRange in standard or
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
 *   boost::for_each(blitzplusplus::basi::fullRange<tmptools::Vector<3,6,1,9,7> >(psi),
 *                   actWithA);
 * }
 * ~~~
 * 
 * \see blitzplusplus::basi::fullRange and the \refBoostConstruct{for_each,range/doc/html/range/reference/algorithms/non_mutating/for_each.html} algorithm of Boost.Range.
 * 
 * For further basic examples of usage cf. `utils/testsuite/BlitzArraySliceIterator.cc` & `utils/testsuite/BlitzArraySliceIteratorTMP.cc`.
 * 
 * \see \ref iteratorimplementation
 * 
 * \todo Implement a default version of Iterator for the case when neither slicing nor transposition is necessary, that is when `V` is equivalent to a range<0,RANK-1>.
 * This will require further compile-time implementation selection.
 *
 * \todo Refine the iterator category according to the \refBoost{New-Style Iterator concepts,iterator/doc/index.html#new-style-iterators}.
 * The point is that a multi-array is not a container of slices, so Iterator is definitely not a standard iterator. It seems rather like a proxy iterator.
 *
 */
template<int RANK, typename V, bool IS_CONST>
class Iterator 
  : public ttd::ForwardIteratorHelper<Iterator<RANK,V,IS_CONST>,V,IS_CONST>, // The inheritance has to be public here because of types necessary to inherit
    public BASE_class,
    private ConsistencyChecker<RANK,V>
{
public:
  typedef typename BASE_class Base;

#undef  BASE_class

  typedef typename Base::CcCA CcCA; ///< ttd::ConditionalConstCArray

  typedef boost::iterator_range<Iterator> Range; ///< Boost.Range-compliant range

  Iterator& operator++() {Base::increment(); return *this;} ///< For the ForwardIterator concept

  /// Can be initialized either to the beginning or the end of the sequence of dummy-index combinations
  /** \tparam IS_END governs the end-ness */
  template<bool IS_END>
  Iterator(CcCA& array, boost::mpl::bool_<IS_END> isEnd) : Base(array,isEnd) {}

};


#define NS_NAME basi
#define RETURN_type1(IS_CONST) Iterator<Rank<A>::value,V_S,IS_CONST>
#define ADDITIONAL_PARAMETER
#define ADDITIONAL_ARGUMENT

#include "details_BlitzArraySliceIteratorReentrant.h"

namespace details {

///////
//
// Base
//
///////

/// Generic worker base for Iterator
template<int RANK, typename V, bool IS_CONST>
class Base : public Indexer<RANK,V> ///< This inheritance is only to facilitate use in LazyDensityOperatorSliceIterator, and similar situations. The same does not appear necessary for BaseSpecial
{
public:
  typedef                      CArray<RANK>            CA   ;
  typedef ttd::ConditionalConstCArray<RANK,IS_CONST> CcCA   ;
  typedef ttd::ResCArray<V>                            CARes;
  typedef ttd::ConditionalConstResCArray<V,IS_CONST> CcCARes;

  typedef ttd::VecIdxTiny<RANK,V> VecIdxTiny;

  static const int RANKIDX=RANK-Size<V>::value;

  typedef cpputils::MultiIndexIterator<RANKIDX> Impl;

  Base(CcCA&, Begin); // if it's not the end, it's the beginning
  Base(CcCA&, End  );

  void increment() {++impl_;}

  CcCARes& operator*() const {return Indexer<RANK,V>::index(array_,arrayRes_,*impl_);}
  // This has to return a reference, cf eg. ElementLiovillean::JumpStrategy takes its argument as a non-const reference. Or, it has to take a copy there...

  friend bool operator==(const Base& i1, const Base& i2) {return i1.impl_==i2.impl_ /* && i1.array_==i2.array_ */;}
  // The user has to ensure that the two arrays are actually the same

  const Impl& operator()() const {return impl_;}

private:
  template<bool TAG>
  static Impl ctorHelper(CcCA&);

  mutable CA array_;
  // By value, so that the necessary transposition specified by V can be readily performed. In this case, however, neither stored array can ever be const.

  mutable CARes arrayRes_;

  Impl impl_;

};


//////////////
//
// BaseTrivial
//
//////////////

class OutOfRange : public cpputils::Exception {};


/// used in the case when `Size<V> = RANK` and when `RANK = 1`
template<typename V, bool IS_CONST>
class BaseTrivial
{
public:
  static const int RANK=Size<V>::value;

  typedef                      CArray<RANK>             CA;
  typedef ttd::ConditionalConstCArray<RANK,IS_CONST>  CcCA;

  BaseTrivial(CcCA&, Begin); // if it's not the end, it's the beginning
  BaseTrivial(CcCA&, End  );

  void increment() {if (!isEnd_) isEnd_=true; else throw OutOfRange();}

  CcCA& operator*() const {if (isEnd_) throw OutOfRange(); return array_;}

  friend bool operator==(const BaseTrivial& i1, const BaseTrivial& i2) {return i1.isEnd_==i2.isEnd_;}

protected:
  mutable CA array_;

private:
  bool isEnd_;

};

//////////////
//
// BaseSpecial
//
//////////////

/// used in the case when `Size<V> = RANK`
template<typename V, bool IS_CONST>
class BaseSpecial : public BaseTrivial<V,IS_CONST>
{
public:
  typedef typename BaseTrivial<V,IS_CONST>::CcCA CcCA;
  
  BaseSpecial(CcCA&, Begin); // if it's not the end, it's the beginning
  BaseSpecial(CcCA&, End  );

};

} // details

} // basi



/// Contains data for pre-calculated slices for basi_fast::Iterator
/**
 * This is most useful in situations where the same structure of slices is needed for several arrays of the same structure or several times for the same array.
 * Then, the data necessary for slicing calculated once, can be used for the different arrays.
 * 
 * Another positive is the easy implementation of basi_fast::Iterator, which basically only needs to iterate over the data of the different slices.
 * Depending on the data structure SlicesData::Impl, basi_fast::Iterator can be of different categories. (E.g. bidirectional iterator for a list,
 * but random-access iterator for a vector.)
 * 
 * The drawback is less convenient usage than that of basi::Iterator, since here a separated instance of SlicesData must be stored.
 * 
 * \tparamRANK \tparamV
 * 
 */
template<int RANK, typename V>
class SlicesData
{
public:
  typedef std::list<ptrdiff_t> Impl; ///< Data structure for the sequence of slices

  friend class basi_fast::Iterator<RANK,V, true>;
  friend class basi_fast::Iterator<RANK,V,false>;

  /// Constructor from a reference array
  /** \param array reference array showing the structure for which the slicing is envisaged */
  SlicesData(const CArray<RANK>& array);

private:
  static const Impl ctorHelper(const CArray<RANK>&);

  const Impl firstOffsets_;

  const blitz::TinyVector<int      ,basi::Size<V>::value>  shape_;
  const blitz::TinyVector<ptrdiff_t,basi::Size<V>::value> stride_;

  const blitz::GeneralArrayStorage<basi::Size<V>::value> storage_;

};



/// Contains a “fast” version of BlitzArraySliceIterator
namespace basi_fast {

/// “Fast” version of basi::Iterator relying on a pre-calculated set of slices stored by SlicesData
/**
 * Structure analogous to basi::Iterator.
 * 
 * \tparamRANK \tparamV
 * 
 * \see SlicesData
 * 
 */
template<int RANK, typename V, bool IS_CONST>
class Iterator 
  : public basi::ttd::ForwardIteratorHelper<Iterator<RANK,V,IS_CONST>,V,IS_CONST>
{
public:
  typedef basi::ttd::ConditionalConstCArray<RANK,IS_CONST> CcCA   ;
  typedef basi::ttd::ResCArray<V>                            CARes;
  typedef basi::ttd::ConditionalConstResCArray<V,IS_CONST> CcCARes;

  typedef boost::iterator_range<Iterator> Range;

  Iterator& operator++() {++iter_; return *this;}

  CcCARes& operator*() const 
  {
    arrayRes_.reference(CARes(arrayData_+*iter_,slicesData_.shape_,slicesData_.stride_,blitz::neverDeleteData,slicesData_.storage_));
    return arrayRes_;
  }

  friend bool operator==(const Iterator& i1, const Iterator& i2) {return i1.iter_==i2.iter_;}

  template<bool IS_END>
  Iterator(CcCA& array, const SlicesData<RANK,V>&, boost::mpl::bool_<IS_END> b);

private:
  typename SlicesData<RANK,V>::Impl::const_iterator iter_;

  mutable CARes arrayRes_;

  dcomp*const arrayData_;
  // This can be non_const because the constructor and the dereferencing operator will anyway ensure const-correctness

  const SlicesData<RANK,V>& slicesData_;

};



#define NS_NAME basi_fast
#define RETURN_type1(IS_CONST) Iterator<Rank<A>::value,V_S,IS_CONST>
#define ADDITIONAL_PARAMETER , sd
#define ADDITIONAL_ARGUMENT  , const SlicesData<Rank<A>::value,V_S>& sd

#include "details_BlitzArraySliceIteratorReentrant.h"


} // basi_fast


namespace basi {

/// Iterator to the beginning of the (non-const) sequence \related Iterator
/** \tparam V the retained index positions \tparam A array type taken itself as a parameter to ease template-parameter inference \return non-const Iterator */
template<typename V, typename A>
const Iterator<Rank<A>::value,V,false>
begin(      A& array );

/// ” for the end \related Iterator
template<typename V, typename A>
const Iterator<Rank<A>::value,V,false>
end  (      A& array );

/// Iterator to the beginning of the (const) sequence \related Iterator
/** \tparam V the retained index positions \tparam A array type taken itself as a parameter to ease template-parameter inference \return const Iterator */
template<typename V, typename A>
const Iterator<Rank<A>::value,V,true>
begin(const A& array );

/// “ for the end \related Iterator
template<typename V, typename A>
const Iterator<Rank<A>::value,V,true>
end  (const A& array );

/// \refBoost{Boost.Range,range/doc/html/range/reference/utilities/iterator_range.html}-compliant full range of slice iterators \related Iterator
/**
 * It corresponds to all the possible combinations of dummy indices (\f$\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},\f$…) \ref retainedindexpositionsdefined "here".
 * \tparam V the retained index positions \tparam A array type taken itself as a parameter to ease template-parameter inference 
 */
template<typename V, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,V,false> >
fullRange(      A& array );

/// ” const version \related Iterator
template<typename V, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,V,true> >
fullRange(const A& array );

} // basi

namespace basi_fast {

/// Same as ::begin but here it returns a “fast” Iterator instance \related Iterator
template<typename V, typename A>
const Iterator<Rank<A>::value,V,false>
begin(      A& array , const SlicesData<Rank<A>::value,V>& sd);

/// ” \related Iterator
template<typename V, typename A>
const Iterator<Rank<A>::value,V,true>
begin(const A& array , const SlicesData<Rank<A>::value,V>& sd);

/// Same as ::end but here it returns a “fast” Iterator instance \related Iterator
template<typename V, typename A>
const Iterator<Rank<A>::value,V,false>
end  (      A& array , const SlicesData<Rank<A>::value,V>& sd);

/// ” \related Iterator
template<typename V, typename A>
const Iterator<Rank<A>::value,V,true>
end  (const A& array , const SlicesData<Rank<A>::value,V>& sd);

/// Same as ::fullRange but here it returns an iterator range of “fast” Iterator instances \related Iterator
template<typename V, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,V,true> >
fullRange(const A& array , const SlicesData<Rank<A>::value,V>& sd);

/// ” \related Iterator
template<typename V, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,V,false> >
fullRange(      A& array , const SlicesData<Rank<A>::value,V>& sd);

} // basi_fast


namespace basi {

/** \page iteratorimplementation Notes on the implementation of Iterator
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
 * To store and manipulate the heterogenous collection of slicing indices like `0,a,2,a,4,5,a,a,8,a,10` in Indexer::index, the `vector` class of the 
 * \refBoost{Boost.Fusion,fusion/doc/html/index.html} library is used.
 * 
 * Iterator is implemented in terms of the above two helper classes. Each Iterator invokes a Transposer::transpose at its construction, and – if `RANK` is larger
 * than `Size<V>::value` – an Indexer::index @ the point of its dereferencing when the actual slicing occurs.
 * 
 * Iterator is a forward iterator, implemented with the help of \refBoostConstruct{forward_iterator_helper,utility/operators.htm#iterator} from Boost.Operator.
 * For this to work, we need to define only 3 operations:
 * -# Comparison for equality
 * -# Prefix increment
 * -# Dereferencing
 * A special implementation is needed when the size of the compile-time vector `V` equals `RANK` because in this case actually no slicing takes place,
 * only transposition. For this, as at several other places in the framework, we apply conditional inheritance: Iterator inherits from either of two classes 
 * (details::Base or details::BaseSpecial), the decision being made @ compile time with the help of `boost::mpl::if_c`.
 * 
 * The iteration over dummy indices is implemented with the help of cpputils::MultiIndexIterator.
 * 
 * A metaprogramming example {#metaprogrammingexample}
 * =========================
 * 
 * In the following we analyse a metaprogramming example typical for the framework: how the compile-time vector `0,3,2,6,4,5,1,9,8,7,10`
 * for the self-transposition in Transposer::transpose is prepared.
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
 * - `Size<V>::value` must not be larger than the arity
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


} // basi


} // blitzplusplus



#endif // CPPQEDCORE_UTILS_BLITZARRAYSLICEITERATOR_H_INCLUDED


