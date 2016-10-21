// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED

#include "ParticleCavityFwd.h"

#include "ParsParticleCavity.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode_.h"
#include "Particle_.h"

#include "Exception.h"
#include "SmartPtr.h"

#include <boost/utility.hpp>


namespace particlecavity {


struct UnotVClassSignDiscrepancy : public cpputils::Exception {};



typedef quantumoperator::Tridiagonal<2> Tridiagonal;

typedef structure::Hamiltonian<2> Hamiltonian;
typedef quantumoperator::TridiagonalHamiltonian<2,true> TridiagonalHamiltonian;


typedef TridiagonalHamiltonian::Tridiagonals Tridiagonals;


class Base
  : public structure::Interaction<2>
{
protected:
  Base(mode::Ptr, particle::Ptr, double uNot, double etaeff);

};

} // particlecavity


class ParticleOrthogonalToCavity 
  : public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
public:
  template<typename MODE, typename PUMPED_PART>
  ParticleOrthogonalToCavity(const MODE& mode, const PUMPED_PART& part, const particlecavity::ParsOrthogonal& p)
    : ParticleOrthogonalToCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part),p.uNot) {}

private:
  ParticleOrthogonalToCavity(mode::Ptr, particle::PtrPumped, double uNot);
  ParticleOrthogonalToCavity(mode::Ptr, particle::PtrPumped, double uNot, double etaeff);

};



class ParticleAlongCavity
  : private boost::base_from_member<const ModeFunction>,
    public particlecavity::Base, public particlecavity::Hamiltonian
{
private:
  typedef boost::base_from_member<const ModeFunction> MF_Base;
  struct ThePrivateOne {};

public:
  const ModeFunction& getMF() const {return MF_Base::member;}

  template<typename MODE, typename PART>
  ParticleAlongCavity(const MODE& mode, const PART& part, const particlecavity::ParsAlong& p, double vClass=0)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part), p, vClass, ThePrivateOne()) {}

  template<typename MODE, typename PART>
  ParticleAlongCavity(const MODE& mode, const PART& part, const particlecavity::ParsAlongGenericPump& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part), p.uNot, p.kCav, p.modeCav, p.etaeff) {}

  // The following two describe the case when there is an additional fixed standing wave ALONG the cavity, in which case the particle must be derived from PumpedParticleBase
  // The particle type is not a template parameter to avoid ambiguity with the previous constructors

  template<typename MODE>
  ParticleAlongCavity(const MODE& mode, const PumpedParticleBase& part, const particlecavity::ParsAlong& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part), p, ThePrivateOne()) {}

  template<typename MODE>
    ParticleAlongCavity(const MODE& mode, const PumpedParticleBase& part, const particlecavity::ParsAlongGenericPump& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part), p.uNot, p.kCav, p.modeCav, p.etaeff) {}

  template<typename MODE>
  ParticleAlongCavity(const MODE& mode, particle::PtrPumped part, const particlecavity::ParsAlong& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode), part, p, ThePrivateOne()) {}

  template<typename MODE>
    ParticleAlongCavity(const MODE& mode, particle::PtrPumped part, const particlecavity::ParsAlongGenericPump& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),part, p.uNot, p.kCav, p.modeCav, p.etaeff) {}

private:

  typedef quantumoperator::TridiagonalHamiltonian<2,true> TridiagonalHamiltonian;
  typedef TridiagonalHamiltonian::Tridiagonals            Tridiagonals;

  ParticleAlongCavity(mode::Ptr, particle::Ptr,       const particlecavity::ParsAlong& p, double vClass, const ThePrivateOne&);
  ParticleAlongCavity(mode::Ptr, particle::PtrPumped, const particlecavity::ParsAlong& p, const ThePrivateOne&);

  ParticleAlongCavity(mode::Ptr, particle::Ptr,       double uNot, size_t kCav, ModeFunctionType, double etaeff);
  ParticleAlongCavity(mode::Ptr, particle::PtrPumped, double uNot, size_t kCav, ModeFunctionType, double etaeff);

  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const;

  const bool isSpecialH_;
  mutable TridiagonalHamiltonian tridiagonalH_;

  mutable quantumoperator::Tridiagonal<2> firstH_, firstHT_;
  mutable quantumoperator::Tridiagonal<1> secondH_, secondHT_;

};


#endif // CPPQEDELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED
