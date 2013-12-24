// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Pars.h"
#include "Qbit.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;
using qbit::Pars; using qbit::ParsLossy; using qbit::ParsPumped; using qbit::ParsPumpedLossy;

namespace pythonext {

PARS_GETTER_SETTER(double, Pars, delta)
PARS_GETTER_SETTER(dcomp, Pars, qbitInit)
PARS_GETTER_SETTER(dcomp, ParsPumped, eta)
PARS_GETTER_SETTER(double, ParsLossy, gamma)

void export_ParsQbit()
{
  scope namespaceScope = qbitNameSpace;
  class_<Pars>
    (
      "Pars",
      init<parameters::ParameterTable&, const std::string&>()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(delta)
    .PARS_PROPERTY(qbitInit)
  ;
  class_<ParsLossy,bases<Pars> >
    (
      "ParsLossy",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(gamma)
  ;
  class_<ParsPumped,bases<Pars> >
    (
      "ParsPumped",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(eta)
  ;
  class_<ParsPumpedLossy,bases<ParsPumped,ParsLossy> >
    (
      "ParsPumpedLossy",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
  ;
}

} // pythonext