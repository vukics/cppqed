// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARTICLETWOMODES__H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARTICLETWOMODES__H_INCLUDED

#include "ParticleCavity_.h"

#include "TridiagonalHamiltonian.h"



class ParticleTwoModes : public structure::Interaction<3>, public structure::Hamiltonian<3>
{
private:
  ParticleTwoModes(mode::Ptr, mode::Ptr, particle::Ptr, double uNot0, double uNot1, const ModeFunction&, const ModeFunction&, double);

public:
  ParticleTwoModes(mode::Ptr m0, mode::Ptr m1, particle::Ptr p,
                   const particlecavity::ParsAlong& p0, const particlecavity::ParsAlong& p1, double phi=0)
    : ParticleTwoModes(m0,m1,p,p0.uNot,p1.uNot,ModeFunction(p0.modeCav,p0.kCav),ModeFunction(p1.modeCav,p1.kCav),phi) {}

private:
  typedef quantumoperator::TridiagonalHamiltonian<3,true> TridiagonalHamiltonian;
  typedef TridiagonalHamiltonian::Tridiagonals            Tridiagonals;

  void addContribution_v(double, const quantumdata::StateVectorLow<3>&, quantumdata::StateVectorLow<3>&, double) const override;

  mutable quantumoperator::Tridiagonal<3> firstH_, firstHT_;

  mutable quantumoperator::Tridiagonal<1> secondH_, secondHT_;

  const bool isSpecialH_;
  mutable TridiagonalHamiltonian specialH_;

};

#endif // CPPQEDELEMENTS_INTERACTIONS_PARTICLETWOMODES__H_INCLUDED
