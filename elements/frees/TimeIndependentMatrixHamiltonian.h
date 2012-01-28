// -*- C++ -*-
#ifndef   ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED
#define   ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED

#include "TimeIndependentMatrixHamiltonianFwd.h"

#include "Averaged.h"
#include "Free.h"
#include "Hamiltonian.h"

#include "CMatrix.h"

class TimeIndependentMatrixHamiltonian : public structure::Free, public structure::Hamiltonian<1,structure::NO_TIME>
{
public:
  typedef linalg::CMatrix CMatrix;

  TimeIndependentMatrixHamiltonian(const CMatrix&);
  // For safety, the matrix will be copied in.

private:
  void addContribution(const StateVectorLow&, StateVectorLow&) const; 

  const CMatrix hamiltonianOverI_;

};


template <int RANK, bool IS_TD>
class TimeIndependentMatrixHamiltonianAveraged : public TimeIndependentMatrixHamiltonian, public structure::averaged::Transferring<1,RANK,IS_TD>
{
public:
  typedef structure::Averaged<RANK,IS_TD> Averaged;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  TimeIndependentMatrixHamiltonianAveraged(const CMatrix& matrix, const Averaged& averaged, const LazyDensityOperator& ldo)
    : TimeIndependentMatrixHamiltonian(matrix), structure::averaged::Transferring<1,RANK,IS_TD>(averaged,ldo) {}

};



#endif // ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED
