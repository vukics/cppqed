// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "PythonExtension.h"
#include "Namespaces.h"

#include "ModeFunction.h"
#include "Pars.h"
#include "ParticleCavity_.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

using particlecavity::ParsAlong; using particlecavity::ParsOrthogonal;

namespace pythonext{

PARS_GETTER_SETTER(double, ParsOrthogonal, uNot)

PARS_GETTER_SETTER(size_t, ParsAlong, kCav)
PARS_GETTER_SETTER(ModeFunctionType, ParsAlong, modeCav)


void export_ParsParticleCavity()
{
  scope namespaceScope = particlecavityNameSpace;

  class_<ParsOrthogonal>
    (
      "ParsOrthogonal",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(uNot)
  ;
  class_<ParsAlong, bases<ParsOrthogonal> >
    (
      "ParsAlong",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(kCav)
    .PARS_PROPERTY(modeCav)
  ;
}

}