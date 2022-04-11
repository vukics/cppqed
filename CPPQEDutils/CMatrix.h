// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{defines the typedef linalg::CMatrix and some helpers}
#ifndef CPPQEDCORE_UTILS_CMATRIX_H_INCLUDED
#define CPPQEDCORE_UTILS_CMATRIX_H_INCLUDED

#include "CVector.h"


/// Contains utilities for linear algebra.
namespace linalg {

/// Complex matrix
typedef CArray<2> CMatrix;

/// Applies matrix `m` on `a` and adds the result to `b`
/**
 * \f[b_i+=\sum_j m_{i,j}\,a_j.\f]
 * 
 * \return reference to result `b`
 * 
 */
CVector& apply(const CVector& a, CVector& b, const CMatrix& m);

/// Calculates two times the real part of a matrix *in place*
CMatrix&
calculateTwoTimesRealPartOfSelf(CMatrix&);
// NEEDS_WORK


} // linalg

#endif // CPPQEDCORE_UTILS_CMATRIX_H_INCLUDED
