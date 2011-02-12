// -*- C++ -*-
#ifndef   _COMPLEX_ARRAY_EXTENSIONS_H
#define   _COMPLEX_ARRAY_EXTENSIONS_H

#include "BlitzArray.h"
#include "MathExtensions.h"

#include<boost/mpl/bool.hpp>


namespace blitzplusplus {


inline double sqrAbs(const dcomp& c) {return mathutils::sqrAbs(c);}

BZ_DECLARE_FUNCTION_RET(sqrAbs,double)



template<int TWO_TIMES_RANK>
inline
void
hermitianConjugateSelf(TTD_CARRAY(TWO_TIMES_RANK)&);


template<int TWO_TIMES_RANK>
inline
const TTD_CARRAY(TWO_TIMES_RANK)
hermitianConjugate(const TTD_CARRAY(TWO_TIMES_RANK)&);


template<int RANK1, int RANK2, bool MULT>
const TTD_CARRAY(RANK1+RANK2)
doDirect(const TTD_CARRAY(RANK1)&, const TTD_CARRAY(RANK2)&, boost::mpl::bool_<MULT>);


} // blitzplusplus

#include "impl/ComplexArrayExtensions.tcc"

#endif // _COMPLEX_ARRAY_EXTENSIONS_H
