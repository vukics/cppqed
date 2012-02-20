// -*- C++ -*-
#ifndef _BLITZ_TINY_EXTENSIONS_H
#define _BLITZ_TINY_EXTENSIONS_H

#include "BlitzTinyExtensionsFwd.h"

#ifndef   NDEBUG
#include "Exception.h"
#endif // NDEBUG

#include<blitz/tinyvec2.h>


namespace blitzplusplus {


template<typename T1, typename T2, int RANK1, int RANK2>
blitz::TinyVector<T1,RANK1+RANK2>
concatenateTinies(const blitz::TinyVector<T1,RANK1>&, const blitz::TinyVector<T2,RANK2>&);


#ifndef   NDEBUG
struct HalfCutTinyFishyException : cpputils::Exception {};
#endif // NDEBUG

template<typename T, int TWO_TIMES_RANK> 
blitz::TinyVector<T,TWO_TIMES_RANK/2> 
halfCutTiny(const blitz::TinyVector<T,TWO_TIMES_RANK>&);


} // blitzplusplus


#include "impl/BlitzTinyExtensions.tcc"


#endif // _BLITZ_TINY_EXTENSIONS_H
