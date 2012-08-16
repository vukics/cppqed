// -*- C++ -*-
#ifndef   UTILS_INCLUDE_BLITZARRAYEXTENSIONS_H_INCLUDED
#define   UTILS_INCLUDE_BLITZARRAYEXTENSIONS_H_INCLUDED

#include "Exception.h"

#include <blitz/array.h>

#include <boost/math/special_functions/fpclassify.hpp>


namespace blitzplusplus {

inline bool isfinite(double d) {return boost::math::isfinite(d);}

BZ_DECLARE_FUNCTION_RET(isfinite,bool)


class NonContiguousStorageException : public cpputils::Exception {};


template<typename T, int RANK>
const blitz::Array<T,1>
unaryArray(const blitz::Array<T,RANK>&);


class BinaryArrayOrderingErrorException : public cpputils::Exception {};


template<typename T, int TWO_TIMES_RANK>
const blitz::Array<T,2>
binaryArray(const blitz::Array<T,TWO_TIMES_RANK>&);


} // blitzplusplus


#endif // UTILS_INCLUDE_BLITZARRAYEXTENSIONS_H_INCLUDED
