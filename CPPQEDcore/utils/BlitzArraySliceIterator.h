// -*- C++ -*-
/// \briefFile{Definition of blitzplusplus::basi::Iterator together with its helpers}

/** \page multiarrayconcept The multi-array concept
 * 
 * ### Synopsis: The state vector as a multi-array
 * 
 * First, we introduce basic definitions on the algebra of composite quantum systems, that is, when the state vector of the system
 * is an element of a Hilbert space which is the direct product of elementary Hilbert spaces 
 * (by elementary we mean that it cannot be further decomposed as a direct product of more elementary Hilbert spaces):
 * \f[\HSpace=\bigotimes_i\HSpace_i,\quad\ket\iota\in\HSpace,\quad\ket{\iota_i}\in\HSpace_i,\quad\ket\iota=\bigotimes_i\ket{\iota_i}\equiv\ket{\iota_0,\iota_1,…}\f]
 * The number of elementary Hilbert spaces (the number of quantum numbers of the system) is referred to throughout as the *rank* or *arity*
 * (un<em>ary,</em> bin<em>ary,</em> tern<em>ary,</em> quatern<em>ary,</em> etc.) of the system.
 * 
 * Via an example we define *state-vector slices*:
 * \anchor retainedindexpositionsdefined
 * \f[\ket\Psi\equiv\sum_\iota\Psi_\iota\ket\iota\in\HSpace,\quad\ket{\Psi^{\avr{1,3,6,7,9}}(\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},…)}\equiv\sum_{\iota_1,\iota_3,\iota_6,\iota_7,\iota_9}\Psi_{\iota}\ket{\iota_1,\iota_3,\iota_6,\iota_7,\iota_9}\in\bigotimes_{i=1,3,6,7,9}\HSpace_i\f]
 * A state-vector slice is defined by the *retained index positions* \f$\avr{1,3,6,7,9}\f$, which define the subsystem, and the <em>“dummy” indices</em> 
 * \f$(\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},…)\f$. In situations when slicing occurs in the framework, the set of retained index positions 
 * is an information available at compile time, while the set of dummy indices is an information becoming available only at runtime.
 * 
 * Slicing is fully recursive in that a state-vector slice behaves exactly as a state vector, only with a lower rank. It can even be further sliced.
 * It is in particular true that
 * \f[\braket\iota\Psi=\braket{\iota_1,\iota_3,\iota_6,\iota_7,\iota_9}{\Psi^{\avr{1,3,6,7,9}}(\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},…)}\f]
 * 
 * Via an example we define *canonical operator extensions*:
 * \f[A\equiv\sum_kA_\text{k,3}\otimes A_\text{k,6}\otimes A_\text{k,1}\otimes A_\text{k,9}\otimes A_\text{k,7}\in\Lfrak\lp\HSpace_\text{3}\otimes\HSpace_\text{6}\otimes\HSpace_\text{1}\otimes\HSpace_\text{9}\otimes\HSpace_\text{7}\rp\f]
 * \f[A^{\avr{3,6,1,9,7}}(\HSpace)\equiv\sum_k\lp\idop_0\otimes A_\text{k,1}\otimes\idop_2\otimes A_\text{k,3}\otimes\idop_4\otimes\idop_5\otimes A_\text{k,6}\otimes A_\text{k,7}\otimes\idop_8\otimes A_\text{k,9}\otimes\idop_{10}…\rp\in\Lfrak(\HSpace)\f]
 * When the numbers in the angular brackets are permutations of a sequence of ordinals, this is in fact not even an extension,
 * only a permutation of the underlying elementary Hilbert spaces.
 * 
 * Matrix elements of the operator in extended Hilbert spaces can then be calculated by acting with the (possibly permutated) 
 * original operator on an appropriate vector slice:
 * \f[\bra\iota A^{\avr{3,6,1,9,7}}(\HSpace)\ket\Psi=\bra{\iota_1,\iota_3,\iota_6,\iota_7,\iota_9}A^{\avr{1,2,0,4,3}}\ket{\Psi^{\avr{1,3,6,7,9}}(\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},…)}\f]
 * 
 * 
 * ### The `Array` class of the Blitz++ library
 * 
 * Due to the abovementioned recursiveness, the state vector of a composite quantum system is most conveniently represented as a complex
 * (dcomp) multi-array. For a definition of the multi-array concept cf. the \refBoost{Boost.MultiArray manual,multi_array/doc/reference.html#MultiArray}.
 * 
 * By virtue of its adequacy for numerics and its efficiency, we have chosen the `Array` class from the [Blitz++ library](http://blitz.sourceforge.net) 
 * to represent state vectors on the lowest level in the framework.
 * 
 * In namespace cpputils, our collection of general-purpose modules, we rely on the template alias CArray, while @ higher levels of the framework,
 * we use the more intuitive name quantumdata::Types::StateVectorLow.
 *
 * \see noteonusingblitz
 * 
 */


#ifndef UTILS_INCLUDE_BLITZARRAYSLICEITERATOR_H_INCLUDED
#define UTILS_INCLUDE_BLITZARRAYSLICEITERATOR_H_INCLUDED

#include "BlitzArraySliceIteratorFwd.h"

#include "BlitzArray.h"
#include "MultiIndexIterator.h"
#include "TMP_Tools.h"

#include "Exception.h"

#include <boost/mpl/sort.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/max_element.hpp>

#include <boost/range.hpp>

#include <list>



namespace blitzplusplus {

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
 * Implemented via the static assertion \refBoostConstruct{BOOST_MPL_ASSERT_MSG,mpl/doc/refmanual/assert.html} from Boost.MPL.
 * 
 * - size of `V` must not be larger than `RANK`
 * - `V` must not contain duplicated elements
 * - all elements of `V` must be smaller than `RANK`
 * 
 * \see \ref specifyingsubsystems
 * 
 */
template<int RANK, typename V> struct ConsistencyChecker
{
  BOOST_MPL_ASSERT_MSG( RANK >= Size<V>::value, INDEXER_with_NONPOSITIVE_RANK, () );

  typedef typename boost::mpl::sort<V>::type SortedV;
  BOOST_MPL_ASSERT_MSG( 
                       (
                        boost::mpl::equal<typename boost::mpl::unique<SortedV,boost::is_same<boost::mpl::_1,boost::mpl::_2> >::type,SortedV>::value 
                        ), 
                       BASI_ITERATOR_INCONSISTENT_VECTOR, (V) 
                        );
  
  BOOST_MPL_ASSERT_MSG(
                       (
                        boost::mpl::deref<typename boost::mpl::max_element<V>::type>::type::value < RANK
                        ),
                       BASI_ITERATOR_VECTOR_OUT_of_RANGE, (V)
                       );
};


/// Specializations necessary since boost::mpl::range_c cannot be sorted
template<int RANK, int I1, int I2>  struct ConsistencyChecker<RANK,boost::mpl::range_c<int,I1,I2> >
{ BOOST_MPL_ASSERT_MSG( ( I1>=0 && I2<RANK ), BASI_ITERATOR_VECTOR_OUT_of_RANGE, ( boost::mpl::range_c<int,I1,I2> ) ) ; };

/** \cond */
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


/// Performs the slicing on an array already transposed by Transposer.
template<int RANK, typename V>
class Indexer : public Transposer<RANK,V>
{
public:
  /// Static worker
  /**
   * \param array The array to be sliced.
   * \param resArray The array storing the result.
   * \param idx Set of dummy indices.
   * 
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
  index(CArray<RANK>& array, ttd::ResCArray<V>& resArray, const ttd::VecIdxTiny<RANK,V>& idx)
  {throw TransposerOrIndexerRankTooHighException<RANK>();}

};

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
 * \tparam V compile-time vector holding the *retained index positions* like \f$\avr{3,6,1,9,7}\f$ \ref retainedindexpositionsdefined "above". (Cf. \ref specifyingsubsystems)
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
 * Sticking to the example \ref retainedindexpositionsdefined "above", assume that the function
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
#define RETURN_type1(IS_CONST) Iterator<ArrayRankTraits<A>::value,V_S,IS_CONST>
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
#define RETURN_type1(IS_CONST) Iterator<ArrayRankTraits<A>::value,V_S,IS_CONST>
#define ADDITIONAL_PARAMETER , sd
#define ADDITIONAL_ARGUMENT  , const SlicesData<ArrayRankTraits<A>::value,V_S>& sd

#include "details_BlitzArraySliceIteratorReentrant.h"


} // basi_fast


namespace basi {

/// Iterator to the beginning of the sequence
//@{
template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,true>
begin(const A& array );

template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,false>
begin(      A& array );
//@}

/// Iterator to the end of the sequence
//@{
template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,true>
end (const A& array );

template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,false>
end (      A& array );
//@}

/// \refBoost{Boost.Range,range/doc/html/range/reference/utilities/iterator_range.html}-compliant full range of slice iterators
/** It corresponds to all the possible combinations of dummy indices (\f$\iota_0,\iota_2,\iota_4,\iota_5,\iota_8,\iota_{10},\f$…) \ref retainedindexpositionsdefined "above". */
//@{
template<typename V, typename A>
const boost::iterator_range<Iterator<ArrayRankTraits<A>::value,V,true> >
fullRange(const A& array );

template<typename V, typename A>
const boost::iterator_range<Iterator<ArrayRankTraits<A>::value,V,false> >
fullRange(      A& array );
//@}

} // basi

namespace basi_fast {

/// Same as blitzplusplus::basi::begin, blitzplusplus::basi::end, and blitzplusplus::basi::fullRange, respectively, but here they return “fast” Iterator instances.
//@{
template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,true>
begin(const A& array , const SlicesData<ArrayRankTraits<A>::value,V>& sd);

template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,true>
end (const A& array , const SlicesData<ArrayRankTraits<A>::value,V>& sd);

template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,false>
begin(     A& array , const SlicesData<ArrayRankTraits<A>::value,V>& sd);

template<typename V, typename A>
const Iterator<ArrayRankTraits<A>::value,V,false>
end (      A& array , const SlicesData<ArrayRankTraits<A>::value,V>& sd);

template<typename V, typename A>
const boost::iterator_range<Iterator<ArrayRankTraits<A>::value,V,true> >
fullRange(const A& array , const SlicesData<ArrayRankTraits<A>::value,V>& sd);

template<typename V, typename A>
const boost::iterator_range<Iterator<ArrayRankTraits<A>::value,V,false> >
fullRange(      A& array , const SlicesData<ArrayRankTraits<A>::value,V>& sd);
//@}

} // basi_fast


namespace basi {

/** \page iteratorimplementation Notes on the implementation of Iterator
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
 * ### A metaprogramming example
 * 
 * In the following we analyse a metaprogramming example typical for the framework: how the compile-time vector `0,3,2,6,4,5,1,9,8,7,10`
 * for the self-transposition in Transposer::transpose is prepared.
 * 
 * This is done by the following snippet in `utils/BlitzArraySliceIterator.tcc`:
 * \snippet utils/BlitzArraySliceIterator.tcc A metaprogramming example
 * \par Line 9:
 * We are using the \refBoostConstruct{fold,mpl/doc/refmanual/fold.html} metaalgorithm from Boost.MPL.
 * Here, it iterates over the sequence of ordinals (tmptools::Ordinals) between `0` and `RANK-1`.
 * \par Line 10:
 * The initial state for the fold algorithm is an empty \refBoost{compile-time vector of integers,mpl/doc/refmanual/vector-c.html}
 * and the \refBoost{iterator,mpl/doc/refmanual/begin.html} pointing to the first element of the compile-time vector `V`. 
 * These two are “zipped” into a \refBoost{compile-time pair,mpl/doc/refmanual/pair.html}.
 * At the end, the first element of this pair will hold the result.
 * \par Lines 11-21
 * express the forward operation of the fold algorithm in the form of a \refBoost{compile-time lambda expression,mpl/doc/tutorial/handling-placeholders.html}.
 * After each step, the new state will again be a pair composed of
 * -# (Lines 11-16:) the vector is \refBoost{augmented,mpl/doc/refmanual/push-back.html} either by the current element in `V`
 *    pointed to by the iterator (Line 13), or the current ordinal (Line 14), depending on whether `V` tmptools::numerical_contains this ordinal.
 * -# (Lines 17-20:) the iterator is \refBoost{advanced,mpl/doc/refmanual/next.html} (Line 18) if the same numerical containment
 *    criterion is met, otherwise it is left untouched for the same element to be considered again (Line 19).
 *    \note Note that in template metaprogramming, “calling” tmptools::numerical_contains<V,mpl::_2> twice is no waste because the second time no further template 
 * instantiations are needed, the compiler will use the ones already instantiated the first time.
 * \par Line 28:
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



#endif // UTILS_INCLUDE_BLITZARRAYSLICEITERATOR_H_INCLUDED


