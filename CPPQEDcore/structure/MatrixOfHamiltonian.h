// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_STRUCTURE_MATRIXOFHAMILTONIAN_H_INCLUDED
#define   CPPQEDCORE_STRUCTURE_MATRIXOFHAMILTONIAN_H_INCLUDED

#include "StateVector.tcc"

// NEEDS_WORK the input of hamiltonian should be sharedPointerized here, in which case RANK could be extracted through a traits class

template<typename T> // T should model both Hamiltonian and DimensionsBookkeeper
const linalg::CMatrix
calculateMatrix(const T& hamiltonian, double t=0, double t0=0)
{
  using namespace linalg;

  static const int RANK=T::DIMESIONS_BOOKKEEPER_RANK;
  typedef typename quantumdata::StateVector<RANK> StateVector;

  const int dim=hamiltonian.getTotalDimension();

  StateVector 
    psi   (hamiltonian.getDimensions()), 
    dpsidt(psi);
  // both are now 0

  CMatrix res(dim,dim);
  
  CVector psi1D(psi.vectorView());

  const CVector dpsidt1D(dpsidt.vectorView());

  // Since the hamiltonian is in general non-Hermitian, all matrix elements have to be separately calculated.
  for (int row=0; row<dim; ++row) {
    psi1D(row)=1;
    for (int col=0; col<dim; ++col) {
      if (const auto ha=dynamic_cast<const structure::Hamiltonian<RANK>*>(&hamiltonian)) 
        ha->addContribution(t,psi.getArray(),dpsidt.getArray(),t0);
      res(col,row)=DCOMP_I*dpsidt1D(col); 
      // note that addContribution is meant to calculate Hamiltonian over DCOMP_I, which we correct here
      dpsidt=0;
    }
    psi1D(row)=0;
  }

  return res;

}




#endif // CPPQEDCORE_STRUCTURE_MATRIXOFHAMILTONIAN_H_INCLUDED
