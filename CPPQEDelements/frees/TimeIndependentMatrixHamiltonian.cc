// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "TimeIndependentMatrixHamiltonian.h"


TimeIndependentMatrixHamiltonian::TimeIndependentMatrixHamiltonian(const CMatrix& hamiltonianOverI)
  : structure::Free(hamiltonianOverI.extent(0),{RF{"Largest frequency",max(abs(hamiltonianOverI)),1.}}),
    hamiltonianOverI_(hamiltonianOverI.copy())
{
  if (hamiltonianOverI_.extent(0)!=hamiltonianOverI_.extent(1)) throw std::invalid_argument("Matrix sot square in TimeIndependentMatrixHamiltonian");
  getParsStream()<<"Time-independent matrix Hamiltonian"<<std::endl;
}


void TimeIndependentMatrixHamiltonian::addContribution_v(structure::NoTime, const structure::freesystem::StateVectorLow& psi, structure::freesystem::StateVectorLow& dpsidt) const
{
  linalg::apply(psi,dpsidt,hamiltonianOverI_);
}
