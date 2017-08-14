// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_TIMEINDEPENDENTMATRIXHAMILTONIAN_H_INCLUDED
#define   CPPQEDELEMENTS_FREES_TIMEINDEPENDENTMATRIXHAMILTONIAN_H_INCLUDED

#include "TimeIndependentMatrixHamiltonianFwd.h"

#include "AveragingUtils.tcc"

#include "Averaged.h"
#include "Free.h"
#include "Hamiltonian.h"

#include "CMatrix.h"

class TimeIndependentMatrixHamiltonian : public structure::Free, public structure::HamiltonianTimeDependenceDispatched<1,structure::NO_TIME>
{
public:
  typedef linalg::CMatrix CMatrix;

  TimeIndependentMatrixHamiltonian(const CMatrix&);
  // For safety, the matrix will be copied in.

private:
  void addContribution_v(structure::NoTime, const StateVectorLow&, StateVectorLow&) const; 

  const CMatrix hamiltonianOverI_;

};


template <int RANK>
class TimeIndependentMatrixHamiltonianAveraged : public TimeIndependentMatrixHamiltonian, public averagingUtils::Transferring<1,RANK,false>
{
private:
  typedef averagingUtils::Transferring<1,RANK,false> Base;
  
public:
  typedef typename Base::AveragedToPtr AveragedToPtr;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  TimeIndependentMatrixHamiltonianAveraged(const CMatrix& matrix, AveragedToPtr averaged, const LazyDensityOperator& ldo)
    : TimeIndependentMatrixHamiltonian(matrix), Base(averaged,ldo) {}

};



#endif // CPPQEDELEMENTS_FREES_TIMEINDEPENDENTMATRIXHAMILTONIAN_H_INCLUDED
