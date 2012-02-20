// -*- C++ -*-
#ifndef _BLITZ_ARRAY_SMART_ITERATOR_H
#define _BLITZ_ARRAY_SMART_ITERATOR_H

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

#define TTD_CONDITIONAL_CONST_CARRAY(RANK,CONST) typename tmptools::ConditionalAddConst<TTD_CARRAY(RANK),(CONST)>::type

#define MPL_SIZE(V) (boost::mpl::size<V>::value)

#define TTD_RES_CARRAY(V) TTD_CARRAY(MPL_SIZE(V))

#define TTD_CONDITIONAL_CONST_RES_CARRAY(V,CONST) typename tmptools::ConditionalAddConst<TTD_RES_CARRAY(V),(CONST)>::type

#define TTD_VEC_IDXTINY(RANK,V) TTD_IDXTINY(RANK-MPL_SIZE(V))
// note that Array::lbound and ubound return a TinyVector<int,...>, so that here we have to use int as well.


namespace blitzplusplus {


namespace basi {


namespace details {


struct EmptyBase {};


template<int RANK, typename V> struct ConsistencyChecker
{
  BOOST_MPL_ASSERT_MSG( RANK >= MPL_SIZE(V), INDEXER_with_NONPOSITIVE_RANK, () );

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
template<int RANK, int N>           struct ConsistencyChecker<RANK,tmptools::Ordinals<N>   > : private ConsistencyChecker<RANK,tmptools::Range<N,0> >             {};



template<int RANK, typename V>
TTD_VEC_IDXTINY(RANK,V)
FilterOut(const TTD_IDXTINY(RANK)&);


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
  TTD_CARRAY(RANK)&
  transpose(TTD_CARRAY(RANK)&) {throw TransposerOrIndexerRankTooHighException();}

};


template<int RANK, typename V>
class Indexer : public Transposer<RANK,V>
{
public:
  static
  TTD_RES_CARRAY(V)&
  index(TTD_CARRAY(RANK)&,
	TTD_RES_CARRAY(V)&,
	const TTD_VEC_IDXTINY(RANK,V)&) {throw TransposerOrIndexerRankTooHighException();}

};




//////////////////////////
//
// BlitzArraySliceIterator
//
//////////////////////////

// Think over: is the following solution based on inheritance to solve
// the specialization problem optimal?  Note: this is NOT allowed
// (specialization depending on a template parameter)

// template<typename V, bool CONST> Iterator<MPL_SIZE(V),V,CONST>;


template<typename, bool>
class IteratorBaseSpecial;

template<int, typename, bool>
class IteratorBase;


#define TTD_FORWARD_ITERATOR_HELPER boost::forward_iterator_helper<Iterator<RANK,V,CONST>,TTD_CONDITIONAL_CONST_RES_CARRAY(V,CONST)>

#define BASE_class boost::mpl::if_c<RANK==MPL_SIZE(V),\
				    IteratorBaseSpecial<V,CONST>,\
				    IteratorBase<RANK,V,CONST> >::type

template<int RANK, typename V, bool CONST>
class Iterator 
  : public TTD_FORWARD_ITERATOR_HELPER, // The inheritance has to be public here because of types necessary to inherit
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
  Iterator(CcCA& array, boost::mpl::bool_<IS_END> b) : Base(array,b) {}

};



#define NS_NAME basi
#define RETURN_type1(CONST) Iterator<ArrayRankTraits<A>::value,V_S,CONST>

#include "details/BlitzArraySliceIteratorReentrant.h"

///////
//
// Base
//
///////


template<int RANK, typename V, bool CONST>
class IteratorBase : public Indexer<RANK,V>
// This inheritance is only to facilitate use in LazyDensityOperatorSliceIterator, and similar situations. The same does not appear necessary for IteratorBaseSpecial
{
public:
  typedef TTD_CARRAY(RANK)                            CA   ;
  typedef TTD_CONDITIONAL_CONST_CARRAY(RANK,CONST)  CcCA   ;
  typedef TTD_RES_CARRAY(V)                           CARes;
  typedef TTD_CONDITIONAL_CONST_RES_CARRAY(V,CONST) CcCARes;

  typedef TTD_VEC_IDXTINY(RANK,V) VecIdxTiny;

  static const int RANKIDX=RANK-MPL_SIZE(V);

  typedef cpputils::MultiIndexIterator<RANKIDX> Impl;

  typedef Indexer<RANK,V> Base;

  IteratorBase(CcCA&, boost::mpl::false_); // if it's not the end, it's the beginning
  IteratorBase(CcCA&, boost::mpl:: true_);

  void increment() {++impl_;}

  CcCARes& operator*() const {return Base::index(array_,arrayRes_,*impl_);}
  // This has to return a reference, cf eg
  // ElementLiovillean::JumpStrategy takes its argument as a non-const
  // reference!!! Or, it has to take a copy there...

  friend bool operator==(const IteratorBase& i1, const IteratorBase& i2)
  {return i1.impl_==i2.impl_ /* && i1.array_==i2.array_ */;}
  // The user has to take care that the two arrays are actually the same

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
// BaseSpecial
//
//////////////

class OutOfRange : public cpputils::Exception {};


template<typename V, bool CONST>
class IteratorBaseSpecial
{
public:
  static const int RANK=MPL_SIZE(V);

  typedef TTD_CARRAY(RANK)                            CA;
  typedef TTD_CONDITIONAL_CONST_CARRAY(RANK,CONST)  CcCA;

  IteratorBaseSpecial(CcCA&, boost::mpl::false_); // if it's not the end, it's the beginning
  IteratorBaseSpecial(CcCA&, boost::mpl:: true_);

  void increment() {if (!isEnd_) isEnd_=true; else throw OutOfRange();}

  CcCA& operator*() const {if (isEnd_) throw OutOfRange(); return array_;}

  friend bool operator==(const IteratorBaseSpecial& i1, const IteratorBaseSpecial& i2) {return i1.isEnd_==i2.isEnd_;}

private:
  mutable CA array_;

  bool isEnd_;

};



} // basi




template<int RANK, typename V>
class SlicesData
{
public:
  typedef SlicesData<RANK,V> type;

  typedef TTD_CARRAY(RANK) CArray;

  typedef std::list<ptrdiff_t> Impl;

  friend class basi_fast::Iterator<RANK,V, true>;
  friend class basi_fast::Iterator<RANK,V,false>;

  SlicesData(const CArray&);

private:
  static const Impl ctorHelper(const CArray&);

  static const V v_;

  const Impl firstOffsets_;

  const blitz::TinyVector<int      ,MPL_SIZE(V)>  shape_;
  const blitz::TinyVector<ptrdiff_t,MPL_SIZE(V)> stride_;

  const blitz::GeneralArrayStorage<MPL_SIZE(V)> storage_;

};




namespace basi_fast {


template<int RANK, typename V, bool CONST>
class Iterator 
  : public TTD_FORWARD_ITERATOR_HELPER
{
public:
  typedef TTD_CONDITIONAL_CONST_CARRAY(RANK,CONST)  CcCA   ;
  typedef TTD_RES_CARRAY(V)                           CARes;
  typedef TTD_CONDITIONAL_CONST_RES_CARRAY(V,CONST) CcCARes;

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



#define RETURN_type1(CONST) Iterator<ArrayRankTraits<A>::value,V,CONST>


template<typename A, typename V>
inline
const RETURN_type1(true )
begin(const A& array, const SlicesData<ArrayRankTraits<A>::value,V>& sd)
{return RETURN_type1(true )(array,sd,boost::mpl::false_());}

template<typename A, typename V>
inline
const RETURN_type1(true )
end  (const A& array, const SlicesData<ArrayRankTraits<A>::value,V>& sd)
{return RETURN_type1(true )(array,sd,boost::mpl:: true_());}

template<typename A, typename V>
inline
const RETURN_type1(false)
begin(      A& array, const SlicesData<ArrayRankTraits<A>::value,V>& sd)
{return RETURN_type1(false)(array,sd,boost::mpl::false_());}

template<typename A, typename V>
inline
const RETURN_type1(false)
end  (      A& array, const SlicesData<ArrayRankTraits<A>::value,V>& sd) 
{return RETURN_type1(false)(array,sd,boost::mpl:: true_());}


#define RETURN_type2(CONST) boost::iterator_range<RETURN_type1(CONST)>

template<typename A, typename V>
inline
const RETURN_type2(true ) 
fullRange(const A& array, const SlicesData<ArrayRankTraits<A>::value,V>& sd) 
{return RETURN_type2(true )(blitzplusplus::basi_fast::begin(array,sd),blitzplusplus::basi_fast::end(array,sd));}

template<typename A, typename V>
inline
const RETURN_type2(false) 
fullRange(      A& array, const SlicesData<ArrayRankTraits<A>::value,V>& sd)
{return RETURN_type2(false)(blitzplusplus::basi_fast::begin(array,sd),blitzplusplus::basi_fast::end(array,sd));}


#undef  RETURN_type2
#undef  RETURN_type1



} // basi_fast


} // blitzplusplus


#include "impl/BlitzArraySliceIterator.tcc"

#undef TTD_FORWARD_ITERATOR_HELPER

#undef TTD_VEC_IDXTINY
#undef TTD_CONDITIONAL_CONST_RES_CARRAY
#undef TTD_RES_CARRAY
#undef MPL_SIZE
#undef TTD_CONDITIONAL_CONST_CARRAY

#endif // _BLITZ_ARRAY_SMART_ITERATOR_H


