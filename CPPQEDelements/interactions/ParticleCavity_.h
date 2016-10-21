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


struct UnotEtaeffSignDiscrepancy : public cpputils::Exception {};



typedef quantumoperator::Tridiagonal<2> Tridiagonal;


typedef quantumoperator::TridiagonalHamiltonian<2,true> TridiagonalHamiltonian;


typedef TridiagonalHamiltonian::Tridiagonals Tridiagonals;


class Base
  : public structure::Interaction<2>
{
protected:
  Base(mode::Ptr, particle::Ptr, double uNot, double etaeff);

};

const Tridiagonal interferic(mode::Ptr mode, particle::Ptr particle, double uNotTimesEtaeff, double uNot, const ModeFunction& mf);

class InterferenceBase
  : private boost::base_from_member<const ModeFunction>, 
    public structure::Interaction<2>, public particlecavity::TridiagonalHamiltonian
{
public:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

  const ModeFunction& getMF() const {return MF_Base::member;}

  InterferenceBase(mode::Ptr, particle::Ptr, double u, size_t kCav, ModeFunctionType);

};


} // particlecavity


/*
class Interference
  : public particlecavity::InterferenceBase
{
public:
  typedef particlecavity::InterferenceBase Base;

  Interference()

};
*/

class ParticleOrthogonalToCavity 
  : public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
public:
  template<typename MODE, typename PUMPED_PART>
  ParticleOrthogonalToCavity(const MODE& mode, const PUMPED_PART& part, const particlecavity::ParsOrthogonal& p)
    : ParticleOrthogonalToCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part),p.uNot) {}

private:
  ParticleOrthogonalToCavity(mode::Ptr, particle::PtrPumped, double uNot);

};



class ParticleAlongCavity 
  : private boost::base_from_member<const ModeFunction>,
    public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
private:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

public:
  const ModeFunction& getMF() const {return MF_Base::member;}

  template<typename MODE, typename PART>
  ParticleAlongCavity(const MODE& mode, const PART& part, const particlecavity::ParsAlong& p, double etaeff=0)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part),p.uNot,p.kCav,p.modeCav,etaeff) {}

  // The following two describe the case when there is an additional fixed standing wave ALONG the cavity, in which case PUMPEDPART must be derived from PumpedParticleBase
  // We write two constructors to avoid ambiguity with the previous one

  template<typename MODE>
  ParticleAlongCavity(const MODE& mode, const PumpedParticleBase& part, const particlecavity::ParsAlong& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),cpputils::sharedPointerize(part),p.uNot,p.kCav,p.modeCav,part.getV_Class()) {}

  template<typename MODE>
  ParticleAlongCavity(const MODE& mode, particle::PtrPumped part, const particlecavity::ParsAlong& p)
    : ParticleAlongCavity(cpputils::sharedPointerize(mode),part,p.uNot,p.kCav,p.modeCav,part->getV_Class()) {}

private:

  ParticleAlongCavity(mode::Ptr, particle::Ptr, double uNot, size_t kCav, ModeFunctionType, double etaeff);

};


#endif // CPPQEDELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED
