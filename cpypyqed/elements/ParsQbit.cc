// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Pars.h"
#include "Qbit.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;
using qbit::Pars; using qbit::ParsLossy; using qbit::ParsPumped; using qbit::ParsPumpedLossy; using qbit::ParsLossyPhaseNoise; using qbit::ParsPumpedLossyPhaseNoise;

namespace pythonext {

PARS_GETTER_SETTER(double, Pars, delta)
PARS_GETTER_SETTER(dcomp, Pars, qbitInit)
PARS_GETTER_SETTER(dcomp, ParsPumped, eta)

template <typename BASE=Pars>
double getgamma(const ParsLossy<BASE> *p) {return p->gamma;}

template <typename BASE=Pars>
void setgamma(ParsLossy<BASE> *p, double v) {p->gamma=v;}

template <typename BASE=ParsLossy<>>
double getgamma_parallel(const ParsLossyPhaseNoise<BASE> *p) {return p->gamma_parallel;}

template <typename BASE=ParsLossy<>>
void setgamma_parallel(ParsLossyPhaseNoise<BASE> *p, double v) {p->gamma_parallel=v;}


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
  class_<ParsLossy<>,bases<Pars> >
    (
      "ParsLossy",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .add_property("gamma", &getgamma<>, &setgamma<>)
  ;
  class_<ParsPumped,bases<Pars> >
    (
      "ParsPumped",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(eta)
  ;
  class_<ParsPumpedLossy,bases<ParsPumped> >
    (
      "ParsPumpedLossy",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .add_property("gamma", &getgamma<ParsPumped>, &setgamma<ParsPumped>)
  ;
  class_<ParsLossyPhaseNoise<>,bases<ParsLossy<>> >
    (
      "ParsLossyPhaseNoise",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .add_property("gamma_parallel", &getgamma_parallel<>, &setgamma_parallel<>)
  ;
  class_<ParsPumpedLossyPhaseNoise,bases<ParsPumpedLossy> >
    (
      "ParsPumpedLossyPhaseNoise",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .add_property("gamma_parallel", &getgamma_parallel<ParsPumpedLossy>, &setgamma_parallel<ParsPumpedLossy>)
  ;
}

} // pythonext
