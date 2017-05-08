// Copyright Raimar Sandner 2012–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "ModeFunction.h"
#include "Pars.h"
#include "ParsParticle.h"
#include "ParticleInitialCondition.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

namespace pythonext{

PARS_GETTER_SETTER(double, particle::Pars, omrec)
PARS_GETTER_SETTER(size_t, particle::Pars, fin)
PARS_GETTER(particle::InitialCondition, particle::Pars, init)
PARS_GETTER_SETTER(int, particle::Pars, hoInitn)

PARS_GETTER_SETTER(double, particle::ParsPumped, vClass)
PARS_GETTER_SETTER(size_t, particle::ParsPumped, kPart)
PARS_GETTER_SETTER(ModeFunctionType, particle::ParsPumped, modePart)

void export_ParsParticle()
{
  scope namespaceScope = particleNameSpace;

  class_<particle::Pars>
    (
      "Pars",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(omrec)
    .PARS_PROPERTY(fin)
    .PARS_RO_PROPERTY(init)
    .PARS_PROPERTY(hoInitn)
  ;
  class_<particle::ParsPumped, bases<particle::Pars> >
    (
      "ParsPumped",
      init<parameters::ParameterTable&, optional<const std::string&> >() 
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(vClass)
    .PARS_PROPERTY(kPart)
    .PARS_PROPERTY(modePart)
  ;

}

} // pythonext