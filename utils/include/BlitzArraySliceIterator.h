// -*- C++ -*-
/// \briefFile{definition of blitzplusplus::basi::Iterator together with its helpers}

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
 * (dcomp) multi-array. For a definition of the multi-array concept cf. the
 * [Boost.MultiArray manual](http://www.boost.org/doc/libs/1_53_0/libs/multi_array/doc/reference.html#MultiArray).
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


template <typename V> struct Size : boost::mpl::size<V> {};

template <int RANK, bool IS_CONST> using ConditionalConstCArray=typename tmptools::ConditionalAddConst<CArray<RANK>,IS_CONST>::type;

template <typename V> using ResCArray=CArray<Size<V>::value>;

template <typename V, bool IS_CONST> using ConditionalConstResCArray=ConditionalConstCArray<Size<V>::value,IS_CONST>;

template <int RANK, typename V, bool IS_CONST> using ForwardIteratorHelper=boost::forward_iterator_helper<Iterator<RANK,V,IS_CONST>,ConditionalConstResCArray<V,IS_CONST> >;

template <int RANK, typename V> using VecIdxTiny=IdxTiny<RANK-Size<V>::value>; // note that Array::lbound and ubound return a TinyVector<int,...>, so that here we have to use int as well.


using cpputils::mii::Begin; using cpputils::mii::End;


namespace details {


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



// Specializations: boost::mpl::range_c cannot be sorted
template<int RANK, int I1, int I2>  struct ConsistencyChecker<RANK,boost::mpl::range_c<int,I1,I2> >
{ BOOST_MPL_ASSERT_MSG( ( I1>=0 && I2<RANK ), BASI_ITERATOR_VECTOR_OUT_of_RANGE, ( boost::mpl::range_c<int,I1,I2> ) ) ; };

template<int RANK, int N, int Nbeg> struct ConsistencyChecker<RANK,tmptools::Range<N,Nbeg> > : private ConsistencyChecker<RANK,boost::mpl::range_c<int,Nbeg,Nbeg+N> > {};
template<int RANK, int N>           struct ConsistencyChecker<RANK,tmptools::Ordinals<N>   > : private ConsistencyChecker<RANK,tmptools::  Range  <N,0>             > {};



template<int RANK, typename V>
const VecIdxTiny<RANK,V>
filterOut(const IdxTiny<RANK>&);


} // details


/////////////////////////////////
//
// The essential part of the work
//
/////////////////////////////////

// These are just rudimentary definitions, the actual definitions are as partial specializations (cf. details/IndexerImplementationsSpecializations.h)
// The functions will throw if 

class TransposerOrIndexerRankTooHighException : public cpputils::Exception {};

template<int RANK, typename V>
class Transposer
{
public:
  static 
  CArray<RANK>&
  transpose(CArray<RANK>&) {throw TransposerOrIndexerRankTooHighException();}

};


template<int RANK, typename V>
class Indexer : public Transposer<RANK,V>
{
public:
  static
  ResCArray<V>&
  index(CArray<RANK>&,
        ResCArray<V>&,
        const VecIdxTiny<RANK,V>&) {throw TransposerOrIndexerRankTooHighException();}

};




//////////////////////////
//
// BlitzArraySliceIterator
//
//////////////////////////

// Think over: is the following solution based on inheritance to solve
// the specialization problem optimal?  Note: this is NOT allowed
// (specialization depending on a template parameter)

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
 * \tparam V compile-time vector holding the *retained index positions* like \f$\avr{3,6,1,9,7}\f$ \ref retainedindexpositionsdefined "above". Example models: tmptools::Vector and `mpl::range_c` from 
 *           [Boost.MPL](http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/range-c.html). `mpl::size<V>::value` must not be larger than `RANK`.
 *           `V` must not “contain” negative values, values not smaller than `RANK`, and duplicate values. These are checked for at compile time,
 *           and any violation is signalled by more or less intelligent compiler errors generated with the help of
 *           [Boost.MPL's static assertions](http://www.boost.org/doc/libs/1_53_0/libs/mpl/doc/refmanual/asserts.html).
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
 * This said, it is never really used directly in the framework, but rather through the maker functions below in standard or
 * [Boost.Range algorithms](http://www.boost.org/doc/libs/1_53_0/libs/range/doc/html/range/reference/algorithms.html).
 * 
 * Quite generally, by iterating through all the combinations of indeces *not* belonging to the given subsystem (dummy indeces) and when dereferenced returning the corresponding slice,
 * it can be used to implement the action of operators in extended (and/or permutated) Hilbert spaces.
 *  
 */
template<int RANK, typename V, bool IS_CONST>
class Iterator 
  : public ForwardIteratorHelper<RANK,V,IS_CONST>, // The inheritance has to be public here because of types necessary to inherit
    public BASE_class,
    private details::ConsistencyChecker<RANK,V>
{
public:
  typedef typename BASE_class Base;

#undef  BASE_class

  typedef typename Base::CcCA CcCA;

  typedef boost::iterator_range<Iterator> Range;

  Iterator& operator++() {Base::increment(); return *this;}

  template<bool IS_END>
  Iterator(CcCA& array, boost::mpl::bool_<IS_END> isEnd) : Base(array,isEnd) {}

};



#define NS_NAME basi
#define RETURN_type1(IS_CONST) Iterator<ArrayRankTraits<A>::value,V_S,IS_CONST>
#define ADDITIONAL_PARAMETER
#define ADDITIONAL_ARGUMENT

#include "details/BlitzArraySliceIteratorReentrant.h"

namespace details {

///////
//
// Base
//
///////


template<int RANK, typename V, bool IS_CONST>
class Base : public Indexer<RANK,V>
// This inheritance is only to facilitate use in LazyDensityOperatorSliceIterator, and similar situations. The same does not appear necessary for BaseSpecial
{
public:
  typedef CArray<RANK>                         CA   ;
  typedef ConditionalConstCArray<RANK,IS_CONST> CcCA   ;
  typedef ResCArray<V>                         CARes;
  typedef ConditionalConstResCArray<V,IS_CONST> CcCARes;

  typedef VecIdxTiny<RANK,V> VecIdxTiny;

  static const int RANKIDX=RANK-Size<V>::value;

  typedef cpputils::MultiIndexIterator<RANKIDX> Impl;

  Base(CcCA&, Begin); // if it's not the end, it's the beginning
  Base(CcCA&, End  );

  void increment() {++impl_;}

  CcCARes& operator*() const {return Indexer<RANK,V>::index(array_,arrayRes_,*impl_);}
  // This has to return a reference, cf eg
  // ElementLiovillean::JumpStrategy takes its argument as a non-const
  // reference!!! Or, it has to take a copy there...

  friend bool operator==(const Base& i1, const Base& i2) {return i1.impl_==i2.impl_ /* && i1.array_==i2.array_ */;}
  // The user has to ensure that the two arrays are actually the same

  const Impl& operator()() const {return impl_;}

private:
  template<bool TAG>
  static Impl ctorHelper(CcCA&);

  mutable CA array_;
  // By value, so that the necessary transposition specified by V can
  // be readily performed.
  // In this case, however, neither stored array can ever be const.

  mutable CARes arrayRes_;

  Impl impl_;

};


//////////////
//
// BaseTrivial
//
//////////////

class OutOfRange : public cpputils::Exception {};


template<typename V, bool IS_CONST>
class BaseTrivial
{
public:
  static const int RANK=Size<V>::value;

  typedef CArray<RANK>                          CA;
  typedef ConditionalConstCArray<RANK,IS_CONST>  CcCA;

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




template<int RANK, typename V>
class SlicesData
{
public:
  typedef SlicesData<RANK,V> type;

  typedef CArray<RANK> CArray;

  typedef std::list<ptrdiff_t> Impl;

  friend class basi_fast::Iterator<RANK,V, true>;
  friend class basi_fast::Iterator<RANK,V,false>;

  SlicesData(const CArray&);

private:
  static const Impl ctorHelper(const CArray&);

  const Impl firstOffsets_;

  const blitz::TinyVector<int      ,basi::Size<V>::value>  shape_;
  const blitz::TinyVector<ptrdiff_t,basi::Size<V>::value> stride_;

  const blitz::GeneralArrayStorage<basi::Size<V>::value> storage_;

};



/// Contains a “fast” version of BlitzArraySliceIterator
namespace basi_fast {


template<int RANK, typename V, bool IS_CONST>
class Iterator 
  : public basi::ForwardIteratorHelper<RANK,V,IS_CONST>
{
public:
  typedef basi::ConditionalConstCArray<RANK,IS_CONST> CcCA   ;
  typedef basi::ResCArray<V>                            CARes;
  typedef basi::ConditionalConstResCArray<V,IS_CONST> CcCARes;

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

#include "details/BlitzArraySliceIteratorReentrant.h"


} // basi_fast


} // blitzplusplus


#endif // UTILS_INCLUDE_BLITZARRAYSLICEITERATOR_H_INCLUDED


