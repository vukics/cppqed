// -*- C++ -*-
#ifndef UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
#define UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED

#include "VectorFromMatrixSliceIteratorFwd.h"

#include "BlitzArraySliceIterator.h"


namespace blitzplusplus {


namespace vfmsi {


////////////////////////////////
//
// VectorFromMatrixSliceIterator
//
////////////////////////////////

struct Left  : boost::mpl::false_ {};
struct Right : boost::mpl:: true_ {};


template<int RANK, typename S>
struct LeftRight 
  : tmptools::Range<RANK,
                    boost::mpl::if_<S,
                                    boost::mpl::int_<RANK>,
                                    boost::mpl::int_<0>
                                    >::type::value> {};
// Meaning that
// V: 0 1 2 3 ... RANK-1       for left and
// V: RANK RANK+1 ... 2*RANK-1 for right


template<int RANK, typename S, bool IS_CONST> using Iterator=basi::Iterator<RANK,LeftRight<RANK/2,S>,IS_CONST>;

#define NS_NAME vfmsi
#define RETURN_type1(IS_CONST) Iterator<Rank<A>::value,V_S,IS_CONST>
#define ADDITIONAL_PARAMETER
#define ADDITIONAL_ARGUMENT

#include "details_BlitzArraySliceIteratorReentrant.h"


template<typename S, typename A>
const Iterator<Rank<A>::value,S,true>
begin(const A& array );

template<typename S, typename A>
const Iterator<Rank<A>::value,S,true>
end (const A& array );

template<typename S, typename A>
const Iterator<Rank<A>::value,S,false>
begin(     A& array );

template<typename S, typename A>
const Iterator<Rank<A>::value,S,false>
end (      A& array );

template<typename S, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,S,true> >
fullRange(const A& array );

template<typename S, typename A>
const boost::iterator_range<Iterator<Rank<A>::value,S,false> >
fullRange(      A& array );

} // vfmsi


namespace basi {

template<int RANK, typename S> struct ConsistencyChecker<RANK,blitzplusplus::vfmsi::LeftRight<RANK/2,S> > {};

} // basi


} // blitzplusplus


#endif // UTILS_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
