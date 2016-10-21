// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "CMatrix.h"


namespace linalg {

CVector& apply(const CVector& a, CVector& b, const CMatrix& m)
{
  using namespace blitz::tensor;
  return b+=sum(m(i,j)*a(j),j);
}


CMatrix&
calculateTwoTimesRealPartOfSelf(CMatrix& matrix)
{
  long dim=matrix.extent(0);

  for (int i=0; i<dim; i++) {
    // Diagonal
    matrix(i,i)=2.*real(matrix(i,i));
    for (int j=i+1; j<dim; j++) 
      // Off-diagonal
      matrix(j,i)=conj(matrix(i,j)=matrix(i,j)+conj(matrix(j,i)));
 
  }
  
  return matrix;

}

} // linalg
