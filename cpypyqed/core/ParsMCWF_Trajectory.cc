// Copyright Raimar Sandner 2012–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "Core.h"

#include "Pars.h"
#include "ParsMCWF_Trajectory.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

using trajectory::ParsStochastic;
using std::string;

namespace pythonext {

PARS_GETTER_SETTER(double,   quantumtrajectory::mcwf::Pars, dpLimit)
PARS_GETTER_SETTER(double,   quantumtrajectory::mcwf::Pars, overshootTolerance)
PARS_GETTER_SETTER(int,      quantumtrajectory::mcwf::Pars, logLevel)

void export_20_ParsMCWF_Trajectory()
{
  scope namespaceScope = quantumtrajectoryNameSpace;
  class_<quantumtrajectory::mcwf::Pars, bases<ParsStochastic> >
    (
      "ParsMCWF",
      "Wrapper of :core:`quantumtrajectory::mcwf::Pars`",
      init<parameters::ParameterTable&, optional<const std::string&> >()
    )
    .PARS_PROPERTY(dpLimit)
    .PARS_PROPERTY(overshootTolerance)
    .PARS_PROPERTY(logLevel)
  ;
}

} // pythonext