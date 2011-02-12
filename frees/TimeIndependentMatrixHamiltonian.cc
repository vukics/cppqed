#include "TimeIndependentMatrixHamiltonian.h"

#include<boost/assign/list_of.hpp>


using namespace boost;
using namespace assign;


class MatrixNotSquareException {};


TimeIndependentMatrixHamiltonian::TimeIndependentMatrixHamiltonian(const CMatrix& hamiltonianOverI)
  : structure::Free(hamiltonianOverI.extent(0),tuple_list_of("Largest frequency",max(abs(hamiltonianOverI)),1.)),
    hamiltonianOverI_(hamiltonianOverI.copy())
{
  if (hamiltonianOverI_.extent(0)!=hamiltonianOverI_.extent(1)) throw MatrixNotSquareException();
  getParsStream()<<"# Time-independent matrix Hamiltonian"<<std::endl;
}


void TimeIndependentMatrixHamiltonian::addContribution(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  linalg::apply(psi,dpsidt,hamiltonianOverI_);
}
