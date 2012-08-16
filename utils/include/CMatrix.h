// -*- C++ -*-
#ifndef UTILS_INCLUDE_CMATRIX_H_INCLUDED
#define UTILS_INCLUDE_CMATRIX_H_INCLUDED

#include "CVector.h"


namespace linalg {


typedef TTD_CARRAY(2) CMatrix;


CVector& apply(const CVector&, CVector&, const CMatrix&);

void
calculateTwoTimesRealPartOfSelf(CMatrix&);
// NEEDS_WORK


} // linalg

#endif // UTILS_INCLUDE_CMATRIX_H_INCLUDED
