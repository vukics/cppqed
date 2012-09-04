// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED
#define   ELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED

#include "ParticleCavityFwd.h"

#include "ParsParticleCavity.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode_.h"
#include "Particle_.h"

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
  Base(mode::SmartPtr, particle::SmartPtr, double uNot, double etaeff);

};



class InterferenceBase
  : private boost::base_from_member<const ModeFunction>, 
    public structure::Interaction<2>, public particlecavity::TridiagonalHamiltonian
{
public:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

  const ModeFunction& getMF() const {return MF_Base::member;}

  InterferenceBase(mode::SmartPtr, particle::SmartPtr, double u, size_t kCav, ModeFunctionType);

};


class POC_Base 
  : public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
public:
  POC_Base(mode::SmartPtr, particle::SmartPtrPumped, double uNot);

};


class PAC_Base 
  : private boost::base_from_member<const ModeFunction>,
    public particlecavity::Base, public particlecavity::TridiagonalHamiltonian
{
public:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

  const ModeFunction& getMF() const {return MF_Base::member;}

  PAC_Base(mode::SmartPtr, particle::SmartPtr, double uNot, size_t kCav, ModeFunctionType, double etaeff);

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

  ParticleOrthogonalToCavity(mode::SmartPtr mode, particle::SmartPtrPumped part, const particlecavity::ParsOrthogonal& p)
    : Base(mode,part,p.uNot) {}

};



class ParticleAlongCavity 
  : public particlecavity::PAC_Base
{
public:
  typedef particlecavity::PAC_Base Base;

  ParticleAlongCavity(mode::SmartPtr mode, particle::SmartPtr part, const particlecavity::ParsAlong& p, double etaeff=0)
    : Base(mode,part,p.uNot,p.kCav,p.modeCav,etaeff) {}

  // The following two describe the case when there is an additional fixed standing wave ALONG the cavity
  ParticleAlongCavity(mode::SmartPtr mode, particle::SmartPtrPumped part, const particlecavity::ParsAlong& p)
    : Base(mode,part,p.uNot,p.kCav,p.modeCav,part->getV_Class()) {}


};


#endif // ELEMENTS_INTERACTIONS_PARTICLECAVITY__H_INCLUDED
