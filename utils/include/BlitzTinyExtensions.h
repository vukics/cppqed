// -*- C++ -*-
#ifndef UTILS_INCLUDE_BLITZTINYEXTENSIONS_H_INCLUDED
#define UTILS_INCLUDE_BLITZTINYEXTENSIONS_H_INCLUDED

#include "BlitzTinyExtensionsFwd.h"

#ifndef   NDEBUG
#include "Exception.h"
#endif // NDEBUG

#include <blitz/tinyvec2.h>


namespace blitzplusplus {


template<typename T1, typename T2, int RANK1, int RANK2>
blitz::TinyVector<T1,RANK1+RANK2>
concatenateTinies(const blitz::TinyVector<T1,RANK1>&, const blitz::TinyVector<T2,RANK2>&);


#ifndef   NDEBUG
struct HalfCutTinyException : cpputils::Exception {};
#endif // NDEBUG

template<typename T, int TWO_TIMES_RANK> 
blitz::TinyVector<T,TWO_TIMES_RANK/2> 
halfCutTiny(const blitz::TinyVector<T,TWO_TIMES_RANK>&);


} // blitzplusplus


#endif // UTILS_INCLUDE_BLITZTINYEXTENSIONS_H_INCLUDED
