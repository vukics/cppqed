// -*- C++ -*-
#ifndef _VECTOR_FROM_MATRIX_SMART_ITERATOR_H
#define _VECTOR_FROM_MATRIX_SMART_ITERATOR_H

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
  : tmptools::RangeMF<RANK,
		      boost::mpl::if_<S,
				      boost::mpl::int_<RANK>,
				      boost::mpl::int_<0>
				      >::type::value> {};
// Meaning that
// V: 0 1 2 3 ... RANK-1       for left and
// V: RANK RANK+1 ... 2*RANK-1 for right



#define TTD_VFMSI(RANK,S,CONST) basi::Iterator<RANK,typename LeftRight<RANK/2,S>::type,CONST>

#define NS_NAME vfmsi
#define ADD_PARAMETER
#define ADD_parameter
#define RETURN_type1(CONST) TTD_VFMSI(A::_bz_rank,V,CONST)

#include "details/BlitzArraySliceIteratorReentrant.h"

#undef  TTD_VFMSI

} // vfmsi


namespace vfmsi_fast {

#define TTD_VFMSI(RANK,S,CONST) basi_fast::Iterator<RANK,typename vfmsi::LeftRight<RANK/2,S>::type,CONST>

#define NS_NAME vfmsi_fast
#define ADD_PARAMETER const SlicesData<A::_bz_rank,typename vfmsi::LeftRight<A::_bz_rank/2,V>::type>&
#define ADD_parameter slicesData,
#define RETURN_type1(CONST) TTD_VFMSI(A::_bz_rank,V,CONST)

#include "details/BlitzArraySliceIteratorReentrant.h"

#undef  TTD_VFMSI

} // vfmsi_fast



} // blitzplusplus


#endif // _VECTOR_FROM_MATRIX_SMART_ITERATOR_H
