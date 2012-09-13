// -*- C++ -*-
#ifndef STRUCTURE_HAMILTONIAN_H_INCLUDED
#define STRUCTURE_HAMILTONIAN_H_INCLUDED

#include "HamiltonianFwd.h"

#include "Types.h"


namespace structure {


template<int RANK>
class Hamiltonian<RANK,TWO_TIME> : private quantumdata::Types<RANK>
{
public:
  typedef boost::shared_ptr<const Hamiltonian> Ptr;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  static  void addContribution(double t,
			       const StateVectorLow& psi, StateVectorLow& dpsidt,
			       double tIntPic0,
			       Ptr hamiltonian,
			       StaticTag=theStaticOne) {if (hamiltonian) hamiltonian->addContribution(t,psi,dpsidt,tIntPic0);} 
  // The (in general, non-Hermitian) Hamiltonian part of evolution.
  // dpsidt+=H*psi/DCOMP_I --- Two things to pay attention to: (1) the
  // function calculates the effect of H/DCOMP_I and not merely H (2)
  // the contribution is ADDED to dpsidt instead of REPLACING.

  virtual ~Hamiltonian() {}

  virtual void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const = 0; 

};


template<int RANK>
class Hamiltonian<RANK,ONE_TIME> : public Hamiltonian<RANK>
{
public:
  typedef Hamiltonian<RANK> Base;

  typedef typename Base::StateVectorLow StateVectorLow;

private:
  void addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
  {
    addContribution(t-tIntPic0,psi,dpsidt);
  }

  virtual void addContribution(double, const StateVectorLow&, StateVectorLow&) const = 0;
};


template<int RANK>
class Hamiltonian<RANK,NO_TIME> : public Hamiltonian<RANK>
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

#endif // STRUCTURE_HAMILTONIAN_H_INCLUDED
