// -*- C++ -*-
#ifndef   STRUCTURE_MATRIX_OF_HAMILTONIAN_INCLUDED
#define   STRUCTURE_MATRIX_OF_HAMILTONIAN_INCLUDED

#include "StateVector.h"



template<typename T> // T should model both Hamiltonian and DimensionsBookkeeper
const linalg::CMatrix
calculateMatrix(const T& hamiltonian, double t=0, double tIntPic0=0)
{
  using namespace linalg;

  static const int RANK=T::N_RANK;
  typedef typename quantumdata::StateVector<RANK> StateVector;

  const int dim=hamiltonian.getTotalDimension();

  StateVector 
    psi   (hamiltonian.getDimensions()), 
    dpsidt(psi);
  // both are now 0

  CMatrix res(dim,dim);
  
  CVector psi1D(psi.vectorView());

  const CVector dpsidt1D(dpsidt.vectorView());

  // Since the hamiltonian is in general non-Hermitian, all matrix
  // elements have to be separately calculated.
  for (int row=0; row<dim; ++row) {
    psi1D(row)=1;
    for (int col=0; col<dim; ++col) {
      structure::Hamiltonian<RANK>::addContribution(t,psi(),dpsidt(),tIntPic0,&hamiltonian);
      res(col,row)=DCOMP_I*dpsidt1D(col); 
      // note that addContribution is meant to calculate Hamiltonian over DCOMP_I, which we correct here
      dpsidt=0;
    }
    psi1D(row)=0;
  }

  return res;

}




#endif // STRUCTURE_MATRIX_OF_HAMILTONIAN_INCLUDED
