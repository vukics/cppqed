// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_IMPL_MLJC_TCC_INCLUDED
#define   ELEMENTS_INTERACTIONS_IMPL_MLJC_TCC_INCLUDED

#include "MLJC.h"

#include "impl/MultiLevel.tcc"
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
  ModeDynamics(const MultiLevelBase<NL>& ml, const ModeBase& mode, const dcomp& g)
    : a      (-g*mode::aop(&mode)),
      adagger((g*mode::aop(&mode)).dagger())
  {
    if (dynamic_cast<const multilevel::Exact<NL>*>(&ml)) throw multilevel::MultiLevelExactNotImplementedException();
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
Base<NL,VC>::Base(const MultiLevelBase<NL>* ml, const ModeBase* mode, const VC& gs) 
  : structure::Interaction<2>(Frees(ml,mode),
			      structure::DynamicsBase::RealFreqs(),
			      complexFreqs(gs)), 
    mds_(boost::fusion::transform(gs,CouplingToModeDynamics(*ml,*mode)))
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
void Base<NL,VC>::addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  boost::fusion::for_each(mds_,ElementaryCoupling(t-tIntPic0,psi,dpsidt));
} 


} // mljc

#endif // ELEMENTS_INTERACTIONS_IMPL_MLJC_TCC_INCLUDED
