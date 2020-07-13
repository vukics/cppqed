// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "PythonExtension.h"

#include "ParticleTwoModes.h"

#include "ParsParticleCavity.h"
#include "Pars.h"

using namespace boost::python;

using particlecavity::ParsAlong; using particlecavity::ParsOrthogonal;

namespace pythonext{

void export_ParticleCavity()
{
  class_<ParticleAlongCavity, bases<structure::Interaction<2> >, boost::noncopyable >
    (
      "ParticleAlongCavity",
      init<mode::Ptr, particle::Ptr, const ParsAlong&>()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4> > >()]
    )
    .def(init<mode::Ptr, particle::Ptr, const ParsAlong&, double>()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4> > >()])
    .def(init<mode::Ptr, particle::PtrPumped, const ParsAlong&>()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4> > >()])
  ;
  class_<ParticleOrthogonalToCavity, bases<structure::Interaction<2> >, boost::noncopyable >
    (
      "ParticleOrthogonalToCavity",
      init<mode::Ptr, particle::PtrPumped, const ParsOrthogonal& >()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4> > >()]
    )
  ;
  class_<ParticleTwoModes, bases<structure::Interaction<3> >, boost::noncopyable >
    (
      "ParticleTwoModes",
      init<mode::Ptr, mode::Ptr, particle::Ptr, const ParsAlong&, const ParsAlong& >()
        [with_custodian_and_ward<1,2, with_custodian_and_ward<1,3, with_custodian_and_ward<1,4, with_custodian_and_ward<1,5, with_custodian_and_ward<1,6> > > > >()]
    )
  ;
}

} // pythonext
