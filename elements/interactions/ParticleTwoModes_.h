// -*- C++ -*-
#ifndef   ELEMENT_PARTICLE_TWO_MODES__INCLUDED
#define   ELEMENT_PARTICLE_TWO_MODES__INCLUDED

#include "ParticleTwoModesFwd.h"

#include "ParticleCavity_.h"

#include "TridiagonalHamiltonian.h"

#include "Exception.h"

#include "details/DispatchFreeType.h"


namespace particletwomodes {

class UnotsSignDiscrepancy : public cpputils::Exception {};


class Base : public structure::Interaction<3>, public structure::Hamiltonian<3>
{
public:
  Base(const ModeBase*, const ModeBase*, const ParticleBase*, double uNot0, double uNot1, const ModeFunction&, const ModeFunction&, double);

private:
  void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const; 

  mutable quantumoperator::Tridiagonal<3> firstH_, firstHT_;

  mutable quantumoperator::Tridiagonal<1> secondH_, secondHT_;

};


} // particletwomodes


class ParticleTwoModes : public particletwomodes::Base
{
public:
  typedef particletwomodes::Base Base;

  template<typename F0, typename F1, typename F2>
  ParticleTwoModes(const F0& f0, const F1& f1, const F2& f2, 
		   const particlecavity::ParsAlong& p0, const particlecavity::ParsAlong& p1, double phi=0)
    : Base(dispatchFreeType(f0),dispatchFreeType(f1),dispatchFreeType(f2),p0.uNot,p1.uNot,ModeFunction(p0.modeCav,p0.kCav),ModeFunction(p1.modeCav,p1.kCav),phi) {}
};

#endif // PARTICLE_TWO_MODES__INCLUDED
