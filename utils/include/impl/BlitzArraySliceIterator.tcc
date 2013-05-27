// -*- C++ -*-
#if !BOOST_PP_IS_ITERATING

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



namespace blitzplusplus {

namespace basi {

namespace details {


template<int RANK, typename V>
class FilterOut
{  
public:
  typedef IdxTiny   <RANK>      IdxTiny;
  typedef VecIdxTiny<RANK,V> VecIdxTiny;

  FilterOut(const IdxTiny& from, VecIdxTiny& to) : from_(from), to_(to), curr_(0) {}

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
const VecIdxTiny<RANK,V>
filterOut(const IdxTiny<RANK>& v)
{
  using namespace boost::mpl;
  using namespace tmptools;

  VecIdxTiny<RANK,V> res;

  {
    FilterOut<RANK,V> helper(v,res);
    for_each<filter_view<Ordinals<RANK>,not_<numerical_contains<V,_> > > >(helper);
  }

  return res;

}


///////////////////////
//
// Ctor implementations
//
///////////////////////


template<int RANK, typename V, bool IS_CONST>
template<bool TAG>
typename Base<RANK,V,IS_CONST>::Impl 
Base<RANK,V,IS_CONST>::ctorHelper(CcCA& array)
{
  return Impl(filterOut<RANK,V>(array.lbound()),filterOut<RANK,V>(array.ubound()),boost::mpl::bool_<TAG>());
}


template<int RANK, typename V, bool IS_CONST>
Base<RANK,V,IS_CONST>::Base(CcCA& array, Begin)
  : array_(), arrayRes_(), impl_(ctorHelper<false>(array))
{
  array_.reference(array);
  Base::transpose(array_);
}


template<int RANK, typename V, bool IS_CONST>
Base<RANK,V,IS_CONST>::Base(CcCA& array, End  ) 
  : array_(), arrayRes_(), impl_(ctorHelper< true>(array))
{
  // Transposition is probably not necessary since this is already the
  // end, and is never dereferenced.
}

template<typename V, bool IS_CONST>
BaseTrivial<V,IS_CONST>::BaseTrivial(CcCA& array, Begin)
  : array_(), isEnd_(false)
{
  array_.reference(array);
}

template<typename V, bool IS_CONST>
BaseTrivial<V,IS_CONST>::BaseTrivial(CcCA&      , End  )
  : array_(), isEnd_(true)
{
}

template<typename V, bool IS_CONST>
BaseSpecial<V,IS_CONST>::BaseSpecial(CcCA& array, Begin)
  : BaseTrivial<V,IS_CONST>(array,cpputils::mii::begin)
{
  Transposer<boost::mpl::size<V>::value,V>::transpose(this->array_);
}


template<typename V, bool IS_CONST>
BaseSpecial<V,IS_CONST>::BaseSpecial(CcCA& array, End  ) 
  : BaseTrivial<V,IS_CONST>(array,cpputils::mii::end)
{
}

} // details

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

  typedef VecIdxTiny<RANK,V> VecIdxTiny;

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
    static ptrdiff_t doIt(const basi::ResCArray<V>& slice, const dcomp* dc)
    {
      return slice.data()-dc;
    }
  };

  RETURN_type res;
  boost::transform(basi::fullRange<V>(array),back_inserter(res),boost::bind(&Helper::doIt,_1,array.data()));
  return res;
}

#undef RETURN_type


template<int RANK, typename V>
SlicesData<RANK,V>::SlicesData(const CArray& array)
  : firstOffsets_(ctorHelper(array)),
    shape_  (basi::begin<V>(array)->shape   ()),
    stride_ (basi::begin<V>(array)->stride  ()),
    storage_(basi::begin<V>(array)->ordering()  ,blitz::TinyVector<bool,basi::Size<V>::value>(true))
{
  assert( ( blitz::all(storage_.ascendingFlag()==blitz::TinyVector<bool,basi::Size<V>::value>(true)) ) );
  assert( ( blitz::all(array   .base         ()==blitz::TinyVector<int ,                RANK>(0   )) ) );
}


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
arrayDataDispatcher(const CArray<RANK>& array)
{
  return const_cast<dcomp*>(array.data());
}

template<int RANK>
inline
dcomp*const
arrayDataDispatcher(      CArray<RANK>& array)
{
  return array.data();
}



} // details


template<int RANK, typename V, bool IS_CONST>
template<bool IS_END>
Iterator<RANK,V,IS_CONST>::Iterator(CcCA& array, const SlicesData<RANK,V>& slicesData, boost::mpl::bool_<IS_END> isEnd)
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


#define BOOST_PP_ITERATION_LIMITS (1,BLITZ_ARRAY_LARGEST_RANK)
#define BOOST_PP_FILENAME_1 "impl/BlitzArraySliceIterator.tcc"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


#endif // UTILS_INCLUDE_IMPL_BLITZARRAYSLICEITERATOR_TCC_INCLUDED


#else  // BOOST_PP_IS_ITERATING


#define rank BOOST_PP_ITERATION()

#define INDEXER_print1(z,m,data) boost::mpl::at_c<TM,m>::type::value
#define INDEXER_print2(z,m,data) boost::fusion::at_c<m>(Base::cache_)

namespace blitzplusplus {
namespace basi {

template<typename V> struct Transposer<rank,V>
{
  typedef typename details::TransposerMeta<rank,V>::type TM;

  typedef CArray<rank> Array;

  static Array& transpose(Array& array)
  {
    array.transposeSelf(BOOST_PP_ENUM(rank,INDEXER_print1,~) );
    return array;
  }

};


template<typename V> struct Indexer<rank,V> : Transposer<rank,V>, private details::IndexerBase<rank,V>
{
  typedef details::IndexerBase<rank,V> Base;

  typedef CArray<rank> Array   ;
  typedef ResCArray<V> ArrayRes;

  static ArrayRes& index(Array& array, ArrayRes& arrayRes, const typename Base::VecIdxTiny& idx)
  {
    Base::fillIdxValues(idx);
    arrayRes.reference(array(BOOST_PP_ENUM(rank,INDEXER_print2,~) ) );
    return arrayRes;
  }

};

}
}

#undef INDEXER_print1
#undef INDEXER_print2

#undef rank


#endif // BOOST_PP_IS_ITERATING
