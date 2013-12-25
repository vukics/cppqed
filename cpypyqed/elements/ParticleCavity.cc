// -*- C++ -*-

#include "PythonExtension.h"

#include "ParticleCavity_.h"

#include "ParsParticleCavity.h"
#include "Mode_.h"
#include "Pars.h"
#include "Particle_.h"

using namespace boost::python;

using particlecavity::ParsAlong; using particlecavity::ParsOrthogonal;

namespace pythonext{

void export_ParticleCavity()
{
  class_<ParticleAlongCavity, bases<structure::Interaction<2> >, boost::noncopyable >
    (
      "ParticleAlongCavity",
      init<const ModeBase&, const ParticleBase&, const ParsAlong& >()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4> > >()]
    )
  ;
  class_<ParticleOrthogonalToCavity, bases<structure::Interaction<2> >, boost::noncopyable >
    (
      "ParticleOrthogonalToCavity",
      init<const ModeBase&, const PumpedParticleBase&, const ParsOrthogonal& >()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4> > >()]
    )
  ;
}

} // pythonext