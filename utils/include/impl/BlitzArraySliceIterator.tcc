// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_BLITZARRAYSLICEITERATOR_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_BLITZARRAYSLICEITERATOR_TCC_INCLUDED

#include "BlitzArraySliceIterator.h"

#include "impl/MultiIndexIterator.tcc"
#include "Range.h"

#include <boost/bind.hpp>

#include <boost/mpl/size.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/filter_view.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/push_back.hpp>

#include <boost/fusion/mpl.hpp> // include all, otherwise difficult to find out what is actually needed
#include <boost/fusion/container/vector/convert.hpp>
#include <boost/fusion/algorithm/iteration/fold.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>


#define DEFINE_BLITZ_ARRAY_SLICE_ITERATOR_MACROS
#include "details/BlitzArraySliceIteratorMacros.h"

namespace blitzplusplus {

namespace basi {

namespace details {


template<int RANK, typename V>
class VecIdxTinyImpl 
{  
public:
  typedef TTD_IDXTINY(RANK)          IdxTiny;
  typedef TTD_VEC_IDXTINY(RANK,V) VecIdxTiny;

  VecIdxTinyImpl(const IdxTiny& from, VecIdxTiny& to) : from_(from), to_(to), curr_(0) {}

  template<typename T>
  void operator()(T) 
  {
    to_(curr_++)=from_(T::value);
  }
  
private:
  const IdxTiny& from_;
  VecIdxTiny& to_;
  unsigned curr_;

};


template<int RANK, typename V>
TTD_VEC_IDXTINY(RANK,V)
FilterOut(const TTD_IDXTINY(RANK)& v)
{
  using namespace boost::mpl;
  using namespace tmptools;

  TTD_VEC_IDXTINY(RANK,V) res;

  {
    VecIdxTinyImpl<RANK,V> body(v,res);
    for_each<filter_view<Ordinals<RANK>,not_<numerical_contains<V,_> > > >(body);
  }

  return res;

}


} // details



///////////////////////
//
// Ctor implementations
//
///////////////////////


template<int RANK, typename V, bool CONST>
template<bool TAG>
typename IteratorBase<RANK,V,CONST>::Impl 
IteratorBase<RANK,V,CONST>::ctorHelper(CcCA& array)
{
  using details::FilterOut;
  return Impl(FilterOut<RANK,V>(array.lbound()),FilterOut<RANK,V>(array.ubound()),boost::mpl::bool_<TAG>());
}


template<int RANK, typename V, bool CONST>
IteratorBase<RANK,V,CONST>::IteratorBase(CcCA& array, boost::mpl::false_)
  : array_(), arrayRes_(), impl_(ctorHelper<false>(array))
{
  array_.reference(array);
  Base::transpose(array_);
}


template<int RANK, typename V, bool CONST>
IteratorBase<RANK,V,CONST>::IteratorBase(CcCA& array, boost::mpl:: true_) 
  : array_(), arrayRes_(), impl_(ctorHelper< true>(array))
{
  // Transposition is probably not necessary since this is already the
  // end, and is never dereferenced.
}


template<typename V, bool CONST>
IteratorBaseSpecial<V,CONST>::IteratorBaseSpecial(CcCA& array, boost::mpl::false_)
  : array_(), isEnd_(false)
{
  array_.reference(array);
  Transposer<boost::mpl::size<V>::value,V>::transpose(array_);
}


template<typename V, bool CONST>
IteratorBaseSpecial<V,CONST>::IteratorBaseSpecial(CcCA&, boost::mpl:: true_) 
  : array_(), isEnd_(true)
{
}



//////////////////////////
//
// Indexer Implementations
//
//////////////////////////

namespace details {


template<int RANK, typename V>
struct IdxTypes : boost::mpl::fold<tmptools::Ordinals<RANK>,
				   boost::fusion::vector<>,
				   boost::mpl::push_back<boost::mpl::_1,
							 boost::mpl::if_<tmptools::numerical_contains<V,boost::mpl::_2>,blitz::Range,int>
							 >
				   >
{};


// This class serves as private base for Indexer, cf. details/IndexerImplementationsSpecializations.h
template<int RANK, typename V>
class IndexerBase
{
protected:
  typedef typename IdxTypes<RANK,V>::type Idx;

  typedef TTD_VEC_IDXTINY(RANK,V) VecIdxTiny;

private:
  typedef typename VecIdxTiny::const_iterator CI;

  struct Helper
  {
    typedef CI result_type;

    const CI operator()(CI iter, blitz::Range& t) const
    {
      t=blitz::Range::all(); return iter;
    }

    const CI operator()(CI iter,          int& t) const
    {
      t=*iter++; return iter;
    }

  };

protected:
  static const Idx&
  fillIdxValues(const VecIdxTiny& idx)
  {
    boost::fusion::fold(cache_,idx.begin(),helper_); return cache_;
  }

  static Idx cache_; 
  // In this way the creation of a default-constructed Helper & Idx objects in fillIdxValues can be avoided

private:
  static const Helper helper_;

};

// Definition of members cache_ (the above is only a declaration!):

template<int RANK, typename V>
typename IndexerBase<RANK,V>::Idx IndexerBase<RANK,V>::cache_;

////////////////////////////
//
// Transpose Implementations
//
////////////////////////////


namespace namehider {

using namespace boost::mpl;
namespace mpl=boost::mpl;
using namespace tmptools;

template<int RANK, typename V>
struct Algorithm 
  : fold<Ordinals<RANK>,
	 pair<vector_c<int>,typename boost::mpl::begin<V>::type>,
	 pair<push_back<mpl::first<mpl::_1>,
			if_<numerical_contains<V,mpl::_2>,
			    deref<second<mpl::_1> >,
			    mpl::_2
			    >
			>,
	      if_<numerical_contains<V,mpl::_2>,
		  mpl::next<second<mpl::_1> >,
		  second<mpl::_1>
		  >
	      >
	 >
{};

} // namehider





template<int RANK, typename V>
struct TransposerMeta : boost::mpl::first<typename namehider::Algorithm<RANK,V>::type>
{};


} // details


} // basi



#define RETURN_type typename SlicesData<RANK,V>::Impl

template<int RANK, typename V>
const RETURN_type
SlicesData<RANK,V>::ctorHelper(const CArray& array)
{
  struct Helper {
    static ptrdiff_t doIt(const TTD_CARRAY(MPL_SIZE(V))& slice, const dcomp* dc)
    {
      return slice.data()-dc;
    }
  };

  RETURN_type res;
  boost::transform(basi::fullRange(array,v_),back_inserter(res),boost::bind(&Helper::doIt,_1,array.data()));
  return res;
}

#undef RETURN_type


template<int RANK, typename V>
SlicesData<RANK,V>::SlicesData(const CArray& array)
  : firstOffsets_(ctorHelper(array)),
    shape_  (basi::begin(array,v_)->shape   ()),
    stride_ (basi::begin(array,v_)->stride  ()),
    storage_(basi::begin(array,v_)->ordering() ,blitz::TinyVector<bool,MPL_SIZE(V)>(true))
{
  assert( ( blitz::all(storage_.ascendingFlag()==blitz::TinyVector<bool,MPL_SIZE(V)>(true)) ) );
  assert( ( blitz::all(array   .base         ()==blitz::TinyVector<int ,       RANK>(0   )) ) );
}


template<int Rank, typename V>
const V SlicesData<Rank,V>::v_ = V();

namespace basi_fast {

namespace details {

#define ITER_DISPATCHER(A1,A2)\
inline \
const std::list<ptrdiff_t>::const_iterator \
iterDispatcher(const std::list<ptrdiff_t>& list, boost::mpl::A1) \
{return list.A2();}

ITER_DISPATCHER(false_,begin)
ITER_DISPATCHER( true_,end  )


template<int RANK>
inline
dcomp*const
arrayDataDispatcher(const TTD_CARRAY(RANK)& array)
{
  return const_cast<dcomp*>(array.data());
}

template<int RANK>
inline
dcomp*const
arrayDataDispatcher(      TTD_CARRAY(RANK)& array)
{
  return array.data();
}



} // details


template<int RANK, typename V, bool CONST>
template<bool IS_END>
Iterator<RANK,V,CONST>::Iterator(CcCA& array, const SlicesData<RANK,V>& slicesData, boost::mpl::bool_<IS_END> isEnd)
  : iter_(details::iterDispatcher(slicesData.firstOffsets_,isEnd)),
    arrayRes_(),
    arrayData_(details::arrayDataDispatcher(array)),
    slicesData_(slicesData)
{
#ifndef   NDEBUG
  for (int i=0; i<RANK; ++i) assert( (array.isRankStoredAscending(i)) );
#endif // NDEBUG
  assert( ( blitz::all(array.base()==blitz::TinyVector<int,RANK>(0)) ) );
}




} // basi_fast


} // blitzplusplus


//////////////////////////////////////////
//
// Brute Force Implementations for Indexer
//
//////////////////////////////////////////


#define BOOST_PP_ITERATION_LIMITS (1,11)
#define BOOST_PP_FILENAME_1 "details/IndexerImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


#define UNDEF_BLITZ_ARRAY_SLICE_ITERATOR_MACROS
#include "details/BlitzArraySliceIteratorMacros.h"

#endif // UTILS_INCLUDE_IMPL_BLITZARRAYSLICEITERATOR_TCC_INCLUDED
