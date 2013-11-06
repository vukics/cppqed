// -*- C++ -*-
#ifndef   INTERACTIONS_PARTICLECAVITY_INTERFERENCEWORKAROUND_H_INCLUDED
#define   INTERACTIONS_PARTICLECAVITY_INTERFERENCEWORKAROUND_H_INCLUDED

#include "ParticleCavity_InterferenceWorkaroundFwd.h"

#include "ParsParticleCavity_InterferenceWorkaround.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.h"

#include "Mode.h" // for the UI
#include "Particle.h"

#include "Exception.h"

#include <boost/utility.hpp>


namespace particlecavity_interferenceworkaround {


struct UnotEtaeffSignDiscrepancy : public cpputils::Exception {};



typedef quantumoperator::Tridiagonal<2> Tridiagonal;


typedef structure::TridiagonalHamiltonian<2,true> TridiagonalHamiltonian;


typedef TridiagonalHamiltonian::Tridiagonals Tridiagonals;

class InterferenceBase
  : private boost::base_from_member<const ModeFunction>, 
    public structure::Interaction<2>, public particlecavity_interferenceworkaround::TridiagonalHamiltonian
{
public:
  typedef boost::base_from_member<const ModeFunction> MF_Base;

  const ModeFunction& getMF() const {return MF_Base::member;}

  InterferenceBase(mode::Ptr, particle::Ptr, double u, size_t kCav, ModeFunctionType);

};

} // particlecavity_interferenceworkaround


class Interference
  : public particlecavity_interferenceworkaround::InterferenceBase
{
public:
  typedef particlecavity_interferenceworkaround::InterferenceBase Base;

  Interference(mode::Ptr mode, particle::Ptr part, /*double u, size_t kCav, ModeFunctionType mf,*/ const particlecavity_interferenceworkaround::ParsInterference& p)
    : InterferenceBase(mode, part, p.uInterference, p.kInterference, p.modeInterference) {}

};

#endif // INTERACTIONS_PARTICLECAVITY_INTERFERENCEWORKAROUND_H_INCLUDED
