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


template <int RANK>
class TimeIndependentMatrixHamiltonianAveraged : public TimeIndependentMatrixHamiltonian, public structure::Averaged<1,false>
{
public:
  typedef structure::Averaged<RANK,true> Averaged;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  TimeIndependentMatrixHamiltonianAveraged(const CMatrix& matrix, const Averaged& averaged, const LazyDensityOperator& ldo)
    : TimeIndependentMatrixHamiltonian(matrix), averaged_(averaged), ldo_(ldo) {}

private:
  void process(Averages& averages) const {averaged_.process(averages);}

  void display(const Averages& averages, std::ostream& os, int n) const {averaged_.display(averages,os,n);}

  void displayKey(std::ostream& os, size_t& n) const {averaged_.displayKey(os,n);}

  size_t nAvr() const {return averaged_.nAvr();}

  const Averages average(const quantumdata::LazyDensityOperator<1>&) const {return averaged_.average(0,ldo_);}

  const Averaged& averaged_;

  const LazyDensityOperator& ldo_;

};



#endif // ELEMENTS_TIME_INDEPENDENT_MATRIX_HAMILTONIAN_INCLUDED
