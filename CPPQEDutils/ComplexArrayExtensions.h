// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Helpers for complex `blitz::Array`s, e.g. Hermitian conjugation of multi-matrices}
#ifndef   CPPQEDCORE_UTILS_COMPLEXARRAYEXTENSIONS_H_INCLUDED
#define   CPPQEDCORE_UTILS_COMPLEXARRAYEXTENSIONS_H_INCLUDED

#include "BlitzArray.h"
#include "CMatrix.h"
#include "MathExtensions.h"

#include<boost/mpl/bool.hpp>


namespace blitzplusplus {


inline double sqrAbs(const dcomp& c) {return mathutils::sqrAbs(c);}

BZ_DECLARE_FUNCTION_RET(sqrAbs,double) ;



template<int TWO_TIMES_RANK>
inline
void
hermitianConjugateSelf(CArray<TWO_TIMES_RANK>&);


template<int TWO_TIMES_RANK>
inline
const CArray<TWO_TIMES_RANK>
hermitianConjugate(const CArray<TWO_TIMES_RANK>&);


/// Direct product/sum
/**
 * Returns the direct product (if `IS_MULTIPLICATION`) \f$A_{i,j}=A1_i*A2_j\f$, or direct sum (otherwise)
 * \f$A_{i,j}=A1_i+A2_j\f$ of `array1` and `array2`, with \f$i,\ j\f$ running through all the multi-indices.
 * The implementation is in terms of blitzplusplus::unaryArray views of the arguments.
 * 
 * \tparam IS_MULTIPLICATION governs whether direct product or sum is to be calculated
*/
template<bool IS_MULTIPLICATION, int RANK1, int RANK2>
const CArray<RANK1+RANK2>
doDirect(const CArray<RANK1>&, const CArray<RANK2>&);


} // blitzplusplus


#endif // CPPQEDCORE_UTILS_COMPLEXARRAYEXTENSIONS_H_INCLUDED
