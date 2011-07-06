// -*- C++ -*-
#ifndef   _PARTICLE_CAVITY_H
#define   _PARTICLE_CAVITY_H

#include "ParticleCavityFwd.h"

#include "ParsParticleCavity.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode.h" // for the UI
#include "Particle.h"

#include "Exception.h"

#include <boost/utility.hpp>


namespace particlecavity {


struct UnotEtaeffSignDiscrepancy : public cpputils::Exception {};



typedef quantumoperator::Tridiagonal<2> Tridiagonal;


typedef structure::TridiagonalHamiltonian<2,true> TridiagonalHamiltonian;


typedef TridiagonalHamiltonian::Tridiagonals Tridiagonals;


class Base
  : public structure::Interaction<2>
{
protected:
  Base(const ModeBase*, const ParticleBase*, double uNot, double etaeff);

};



class InterferenceBase
  : private boost::base_from_member<const ModeFunction>, 
    public structure::Interaction<2>, public particlecavity::TridiagonalHamiltonian
{
public:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

  const ModeFunction& getMF() const {return MF_Base::member;}

  InterferenceBase(const ModeBase*, const ParticleBase*, double u, size_t kCav, ModeFunctionType);

};


class POC_Base 
  : public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
public:
  POC_Base(const ModeBase*, const PumpedParticleBase*, double uNot);

};


class PAC_Base 
  : private boost::base_from_member<const ModeFunction>,
    public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
public:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

  const ModeFunction& getMF() const {return MF_Base::member;}

  PAC_Base(const ModeBase*, const ParticleBase*, double uNot, size_t kCav, ModeFunctionType, double etaeff);

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
  : public particlecavity::POC_Base
{
public:
  typedef particlecavity::POC_Base Base;

  ParticleOrthogonalToCavity(const ModeBase& mode, const PumpedParticleBase& part, const particlecavity::ParsOrthogonal& p)
    : Base(&mode,&part,p.uNot) {}

  ParticleOrthogonalToCavity(mode::SmartPtr  mode, particle::SmartPtrPumped  part, const particlecavity::ParsOrthogonal& p)
    : Base(mode.get(),part.get(),p.uNot) {}

};



class ParticleAlongCavity 
  : public particlecavity::PAC_Base
{
public:
  typedef particlecavity::PAC_Base Base;

  ParticleAlongCavity(const ModeBase& mode, const ParticleBase& part, const particlecavity::ParsAlong& p, double etaeff=0)
    : Base(&mode,&part,p.uNot,p.kCav,p.modeCav,etaeff) {}

  ParticleAlongCavity(mode::SmartPtr  mode, particle::SmartPtr  part, const particlecavity::ParsAlong& p, double etaeff=0)
    : Base(mode.get(),part.get(),p.uNot,p.kCav,p.modeCav,etaeff) {}

  // The following two describe the case when there is an additional fixed standing wave ALONG the cavity
  ParticleAlongCavity(const ModeBase& mode, const PumpedParticleBase& part, const particlecavity::ParsAlong& p)
    : Base(&mode,&part,p.uNot,p.kCav,p.modeCav,part.getV_Class()) {}

  ParticleAlongCavity(mode::SmartPtr  mode, particle::SmartPtrPumped  part, const particlecavity::ParsAlong& p)
    : Base(mode.get(),part.get(),p.uNot,p.kCav,p.modeCav,part.get()->getV_Class()) {}


};


#endif // _PARTICLE_CAVITY_H
