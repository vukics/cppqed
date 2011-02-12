// -*- C++ -*-
#ifndef   _BLITZ_ARRAY_EXTENSIONS_H
#define   _BLITZ_ARRAY_EXTENSIONS_H

#include "Exception.h"

#include<blitz/array.h>


namespace blitzplusplus {


class      NonContiguousStorageException : public cpputils::Exception {};


template<typename T, int RANK>
const blitz::Array<T,1>
rankOneArray(const blitz::Array<T,RANK>&);


class RankTwoArrayOrderingErrorException : public cpputils::Exception {};


template<typename T, int TWO_TIMES_RANK>
const blitz::Array<T,2>
rankTwoArray(const blitz::Array<T,TWO_TIMES_RANK>&);


} // blitzplusplus

#include "impl/BlitzArrayExtensions.tcc"


#endif // _BLITZ_ARRAY_EXTENSIONS_H
