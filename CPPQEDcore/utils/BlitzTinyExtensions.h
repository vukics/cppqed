// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Helpers for `blitz::TinyVector`s}
#ifndef CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_H_INCLUDED
#define CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_H_INCLUDED

#include "BlitzTinyExtensionsFwd.h"

#ifndef   NDEBUG
#include "Exception.h"
#endif // NDEBUG

#include <blitz/tinyvec2.h>


namespace blitzplusplus {


/// Concatenates `tiny1` and `tiny2`
/**
 * Implemented with the help of the compile-time–runtime facility \refBoostConstruct{for_each,mpl/doc/refmanual/for-each.html} from Boost.MPL.
 * 
 * \tparam T1 basic type of one of the operand tiny vectors *and the result*
 * \tparam T2 basic type of the other operand tiny vector. Must be convertible to T1.
 * 
 */
template<typename T1, typename T2, int RANK1, int RANK2>
blitz::TinyVector<T1,RANK1+RANK2>
concatenateTinies(const blitz::TinyVector<T1,RANK1>& tiny1, const blitz::TinyVector<T2,RANK2>& tiny2);


#ifndef   NDEBUG
/// Exception class thrown by halfCutTiny
struct HalfCutTinyException : cpputils::Exception {};
#endif // NDEBUG


/// Returns the first half of a tiny vector containing two equal halves
/**
 * In debug mode, the equality of the two halves is checked, and a HalfCutTinyException is thrown @ violation.
 * 
 * \tparam T basic type of the operand tiny vector
 * \tparam TWO_TIMES_RANK length of the ”
 * 
 */
template<typename T, int TWO_TIMES_RANK> 
blitz::TinyVector<T,TWO_TIMES_RANK/2> 
halfCutTiny(const blitz::TinyVector<T,TWO_TIMES_RANK>&);


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_BLITZTINYEXTENSIONS_H_INCLUDED
