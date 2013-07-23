// -*- C++ -*-
#ifndef   ELEMENTS_FREES_TIMEINDEPENDENTMATRIXHAMILTONIAN_H_INCLUDED
#define   ELEMENTS_FREES_TIMEINDEPENDENTMATRIXHAMILTONIAN_H_INCLUDED

#include "TimeIndependentMatrixHamiltonianFwd.h"

#include "AveragingUtils.h"

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
  void addContribution_v(const StateVectorLow&, StateVectorLow&) const; 

  const CMatrix hamiltonianOverI_;

};


template <int RANK, bool IS_TIME_DEPENDENT>
class TimeIndependentMatrixHamiltonianAveraged : public TimeIndependentMatrixHamiltonian, public averagingUtils::Transferring<1,RANK,IS_TIME_DEPENDENT>
{
public:
  typedef typename structure::Averaged<RANK,IS_TIME_DEPENDENT>::Ptr AveragedPtr;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  TimeIndependentMatrixHamiltonianAveraged(const CMatrix& matrix, AveragedPtr averaged, const LazyDensityOperator& ldo)
    : TimeIndependentMatrixHamiltonian(matrix), averagingUtils::Transferring<1,RANK,IS_TIME_DEPENDENT>(averaged,ldo) {}

};



#endif // ELEMENTS_FREES_TIMEINDEPENDENTMATRIXHAMILTONIAN_H_INCLUDED
