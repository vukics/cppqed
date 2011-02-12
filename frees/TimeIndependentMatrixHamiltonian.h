// -*- C++ -*-
#ifndef   ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED
#define   ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED

#include "TimeIndependentMatrixHamiltonianFwd.h"

#include "Free.h"
#include "TimeIndependentHamiltonian.h"

#include "CMatrix.h"

class TimeIndependentMatrixHamiltonian : public structure::Free, public structure::TimeIndependentHamiltonian<1>
{
public:
  typedef linalg::CMatrix CMatrix;

  TimeIndependentMatrixHamiltonian(const CMatrix&);
  // For safety, the matrix will be copied in.

private:
  void addContribution(const StateVectorLow&, StateVectorLow&) const; 

  const CMatrix hamiltonianOverI_;

};

#endif // ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED
