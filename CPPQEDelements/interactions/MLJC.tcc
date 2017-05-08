// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_MLJC_TCC_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_MLJC_TCC_INCLUDED

#include "MLJC.h"

#include "MultiLevel.tcc"
#include "Sigma.h"

using quantumoperator::Sigma;


///////////////////////////////
//
// Helpers to construct & store
//
///////////////////////////////

namespace mljc {


// #define BASE_class boost::tuple<const Frequencies,Tridiagonal,Tridiagonal>

// Unfortunately, cannot use boost::tuple here since it does not support mutable  

template<int NL, typename VC> template<int N1, int N2>
struct Base<NL,VC>::ModeDynamics : public tmptools::pair_c<N1,N2>
{
  ModeDynamics(MultiLevelPtr ml, mode::Ptr mode, const dcomp& g)
    : a      (-g*mode::aop(mode)),
      adagger((g*mode::aop(mode)).dagger())
  {
    if (dynamic_cast<const multilevel::Exact<NL>*>(ml.get())) throw multilevel::MultiLevelExactNotImplementedException();
  }

  mutable Tridiagonal a, adagger;

};

// #undef  BASE_class


//////////////
//
// Constructor
//
//////////////


#define RETURN_type structure::DynamicsBase::ComplexFreqs

template<typename VP>
const RETURN_type
complexFreqs(const VP& etas)
{
  using namespace std;
  RETURN_type res;
  for_each(etas,multilevel::ElementaryComplexFreqs(res,"g"));
  return res;
}

#undef  RETURN_type


template<int NL, typename VC>
Base<NL,VC>::Base(MultiLevelPtr ml, mode::Ptr mode, const VC& gs) 
  : structure::Interaction<2>(Frees(ml,mode),
			      structure::DynamicsBase::RealFreqs(),
			      complexFreqs(gs)), 
    mds_(boost::fusion::transform(gs,CouplingToModeDynamics(ml,mode)))
{
  getParsStream()<<"# Multi-Level Jaynes-Cummings\n";
}


} // mljc



/////////////////////////
//
// Helper for Hamiltonian
//
/////////////////////////

namespace mljc {


template<int NL, typename VC>
class Base<NL,VC>::ElementaryCoupling
{
public:
  ElementaryCoupling(double dt, const StateVectorLow& psi, StateVectorLow& dpsidt)
    : dt_(dt), psi_(psi), dpsidt_(dpsidt) {}

  template<typename MD>
  void operator()(const MD& modeDynamics) const
  {
    Sigma<MD::first,MD::second> sigma;

    Tridiagonal& 
      a      (modeDynamics.a      .propagate(dt_)),
      adagger(modeDynamics.adagger.propagate(dt_));

    // Note that the multiplication with -g and conj(g) has already
    // been taken care of by ModeDynamics above
    
    (sigma         *adagger).apply(psi_,dpsidt_);
    (sigma.dagger()*a      ).apply(psi_,dpsidt_);

  }

private:
  double dt_;

  const StateVectorLow& psi_;
  StateVectorLow& dpsidt_;

};


template<int NL, typename VC>
void Base<NL,VC>::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const
{
  boost::fusion::for_each(mds_,ElementaryCoupling(t-t0,psi,dpsidt));
} 


} // mljc

#endif // CPPQEDELEMENTS_INTERACTIONS_MLJC_TCC_INCLUDED
