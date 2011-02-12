#include "CMatrix.h"

namespace linalg {

CVector& apply(const CVector& psi, CVector& dpsidt, const CMatrix& M)
{
  using namespace blitz::tensor;
  return dpsidt+=sum(M(i,j)*psi(j),j);
}


void
calculateTwoTimesRealPartOfSelf(CMatrix& matrix)
// Calculates two times the real part of matrix IN PLACE
{
  long dim=matrix.extent(0);

  for (int i=0; i<dim; i++) {
    // Diagonal
    matrix(i,i)=2.*real(matrix(i,i));
    for (int j=i+1; j<dim; j++) 
      // Off-diagonal
      matrix(j,i)=conj(matrix(i,j)=matrix(i,j)+conj(matrix(j,i)));
 
  }

}

} // linalg
