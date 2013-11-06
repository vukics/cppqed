// -*- C++ -*-
#ifndef   INTERACTIONS_PARTICLETWOMODES__H_INCLUDED
#define   INTERACTIONS_PARTICLETWOMODES__H_INCLUDED

#include "ParticleTwoModesFwd.h"

#include "ParticleCavity_.h"

#include "TridiagonalHamiltonian.h"

#include "Exception.h"
#include "SmartPtr.h"


namespace particletwomodes {

class UnotsSignDiscrepancy : public cpputils::Exception {};

} // particletwomodes


class ParticleTwoModes : public structure::Interaction<3>, public structure::Hamiltonian<3>
{
private:
  ParticleTwoModes(mode::Ptr, mode::Ptr, particle::Ptr, double uNot0, double uNot1, const ModeFunction&, const ModeFunction&, double);

public:
  template<typename F0, typename F1, typename F2>
  ParticleTwoModes(const F0& f0, const F1& f1, const F2& f2, 
                   const particlecavity::ParsAlong& p0, const particlecavity::ParsAlong& p1, double phi=0)
    : ParticleTwoModes(cpputils::sharedPointerize(f0),cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f2),
                       p0.uNot,p1.uNot,ModeFunction(p0.modeCav,p0.kCav),ModeFunction(p1.modeCav,p1.kCav),phi) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const; 

  mutable quantumoperator::Tridiagonal<3> firstH_, firstHT_;

  mutable quantumoperator::Tridiagonal<1> secondH_, secondHT_;

};

#endif // INTERACTIONS_PARTICLETWOMODES__H_INCLUDED
