// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFile{Declares extensions for creating vector & matrix views of `blitz::Array`s}
#ifndef   CPPQEDCORE_UTILS_BLITZARRAYEXTENSIONS_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZARRAYEXTENSIONS_H_INCLUDED

#include "Exception.h"

#include <blitz/array.h>

#include <boost/math/special_functions/fpclassify.hpp>


/// Comprises our own extensions to Blitz++
namespace blitzplusplus {

inline bool isfinite(double d) {return boost::math::isfinite(d);}

BZ_DECLARE_FUNCTION_RET(isfinite,bool) ;


inline double selectNegative(double d) {return d<0 ? d : 0;}

BZ_DECLARE_FUNCTION_RET(selectNegative,double) ;


/// Exception class thrown by unaryArray & binaryArray
class NonContiguousStorageException : public cpputils::Exception {};


/// Returns a unary view of `array`
/**
 * This is meant to be used only if the underlying storage is contiguous, because while a multi-rank array may be able to represent
 * a view of memory of some more or less intricate structure pattern (e.g. slices), a unary array is not capable of this.
 * In debug mode, violation is detected at runtime via the `array.isStorageContiguous()` member function and an exception of type
 * NonContiguousStorageException is thrown.
 * 
 * \note The returned array does not take part in the reference-counting mechanism of blitz::Array, therefore, it does not own its data, and …
 * 
 * \warning … the returned array has no means of ensuring that the referenced data still exists.
 * 
 * \tparam T the type the array is based on
 * \tparam RANK the arity of the vector space
 * 
 */
template<typename T, int RANK>
const blitz::Array<T,1>
unaryArray(const blitz::Array<T,RANK>&);


/// Exception class thrown by binaryArray
class BinaryArrayOrderingErrorException : public cpputils::Exception {};


/// Returns a binary view of `array`. `TWO_TIMES_RANK` must be an even number
/**
 * Violation is detected @ compile time by tmptools::AssertEvenAndDivideBy2.
 * 
 * The same requirement of contiguity an the same warning applies as for unaryArray, and in addition, further assumptions
 * on the storage order must be made: The storage of the two multi-indices must not be intertwined and must be layed out
 * in the same way, so that e.g. for `RANK=4`, the member function `array.ordering()` should return an octary tiny vector like:
 * 
 *     <1 3 2 0 | 5 7 6 4>
 * 
 * Violation is detected at runtime, and an exception of type BinaryArrayOrderingErrorException is thrown.
 * 
 * \tparam T the type the array is based on
 * \tparam TWO_TIMES_RANK `2*RANK`, `RANK` being the arity of the vector space
 * 
 */
template<typename T, int TWO_TIMES_RANK>
const blitz::Array<T,2>
binaryArray(const blitz::Array<T,TWO_TIMES_RANK>& array);


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_BLITZARRAYEXTENSIONS_H_INCLUDED
