// -*- C++ -*-
#ifndef UTILS_INCLUDE_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED
#define UTILS_INCLUDE_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED

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



#define TTD_VFMSI(RANK,S,CONST) basi::Iterator<RANK,LeftRight<RANK/2,S>,CONST>

#define NS_NAME vfmsi
#define RETURN_type1(CONST) TTD_VFMSI(ArrayRankTraits<A>::value,V_S,CONST)
#define ADDITIONAL_PARAMETER
#define ADDITIONAL_ARGUMENT

#include "details/BlitzArraySliceIteratorReentrant.h"

#undef  TTD_VFMSI

} // vfmsi


namespace basi {

namespace details {

template<int RANK, typename S> struct ConsistencyChecker<RANK,blitzplusplus::vfmsi::LeftRight<RANK/2,S> > {};

} // details

} // basi


} // blitzplusplus


#endif // UTILS_INCLUDE_VECTORFROMMATRIXSLICEITERATOR_H_INCLUDED