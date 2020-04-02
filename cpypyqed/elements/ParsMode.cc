// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Pars.h"
#include "ParsMode.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

namespace pythonext{

PARS_GETTER_SETTER(size_t, mode::Pars, cutoff)
PARS_GETTER_SETTER(size_t, mode::Pars, minitFock)
PARS_GETTER_SETTER(dcomp, mode::Pars, minit)
PARS_GETTER_SETTER(size_t, mode::Pars, displayLevel)
PARS_GETTER_SETTER(double, mode::Pars, delta)

PARS_GETTER_SETTER(dcomp, mode::ParsPumped, eta)

PARS_GETTER_SETTER(double, mode::ParsLossy, kappa)
PARS_GETTER_SETTER(double, mode::ParsLossy, nTh)

void export_ParsMode()
{
  scope namespaceScope = modeNameSpace;
  class_<mode::Pars>
    (
      "Pars",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(cutoff)
    .PARS_PROPERTY(minitFock)
    .PARS_PROPERTY(minit)
    .PARS_PROPERTY(displayLevel)
    .PARS_PROPERTY(delta)
  ;
  class_<mode::ParsPumped, bases<mode::Pars> >
    (
      "ParsPumped",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(eta)
  ;
  class_<mode::ParsLossy, bases<mode::Pars> >
    (
      "ParsLossy",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(kappa)
    .PARS_PROPERTY(nTh)
  ;
  class_<mode::ParsPumpedLossy,bases<mode::ParsPumped,mode::ParsLossy> >
    (
      "ParsPumpedLossy",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
  ;

}

} // pythonext