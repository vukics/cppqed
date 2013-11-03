#include "TimeIndependentMatrixHamiltonian.h"


class MatrixNotSquareException {};


TimeIndependentMatrixHamiltonian::TimeIndependentMatrixHamiltonian(const CMatrix& hamiltonianOverI)
  : structure::Free(hamiltonianOverI.extent(0),FREQS("Largest frequency",max(abs(hamiltonianOverI)),1.)),
    hamiltonianOverI_(hamiltonianOverI.copy())
{
  if (hamiltonianOverI_.extent(0)!=hamiltonianOverI_.extent(1)) throw MatrixNotSquareException();
  getParsStream()<<"# Time-independent matrix Hamiltonian"<<std::endl;
}


void TimeIndependentMatrixHamiltonian::addContribution_v(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  linalg::apply(psi,dpsidt,hamiltonianOverI_);
}
