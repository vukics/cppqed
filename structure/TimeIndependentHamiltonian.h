// -*- C++ -*-
#ifndef _TIME_INDEPENDENT_HAMILTONIAN_H
#define _TIME_INDEPENDENT_HAMILTONIAN_H

#include "TimeIndependentHamiltonianFwd.h"

#include "Hamiltonian.h"


namespace structure {


template<int RANK>
class TimeIndependentHamiltonian : public Hamiltonian<RANK>
{
public:
  typedef Hamiltonian<RANK> Base;

  typedef typename Base::StateVectorLow StateVectorLow;

private:
  void addContribution(double, const StateVectorLow& psi, StateVectorLow& dpsidt, double) const 
  {
    addContribution(psi,dpsidt);
  }

  virtual void addContribution(const StateVectorLow&, StateVectorLow&) const = 0;
};


} // structure


#endif // _TIME_INDEPENDENT_HAMILTONIAN_H
