// -*- C++ -*-
#ifndef _STRUCTURE_HAMILTONIAN_H
#define _STRUCTURE_HAMILTONIAN_H

#include "HamiltonianFwd.h"

#include "Types.h"


namespace structure {


template<int RANK>
class Hamiltonian : public quantumdata::Types<RANK>
{
public:
  typedef quantumdata::Types<RANK> Base;

  typedef typename Base::StateVectorLow StateVectorLow;

  static  void addContribution(double t,
			       const StateVectorLow& psi, StateVectorLow& dpsidt,
			       double tIntPic0,
			       const Hamiltonian* hamiltonian,
			       StaticTag=theStaticOne) {if (hamiltonian) hamiltonian->addContribution(t,psi,dpsidt,tIntPic0);} 
  // The (in general, non-Hermitian) Hamiltonian part of evolution.
  // dpsidt+=H*psi/DCOMP_I --- Two things to pay attention to: (1) the
  // function calculates the effect of H/DCOMP_I and not merely H (2)
  // the contribution is ADDED to dpsidt instead of REPLACING.

  virtual ~Hamiltonian() {}

  virtual void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const = 0; 

};


} // structure

#endif // _STRUCTURE_HAMILTONIAN_H