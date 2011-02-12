// -*- C++ -*-
#ifndef _C_MATRIX_H
#define _C_MATRIX_H

#include "CVector.h"


namespace linalg {


typedef TTD_CARRAY(2) CMatrix;


CVector& apply(const CVector&, CVector&, const CMatrix&);

void
calculateTwoTimesRealPartOfSelf(CMatrix&);
// NEEDS_WORK


} // linalg

#endif // _C_MATRIX_H
