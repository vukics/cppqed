// Copyright Raimar Sandner 2012–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#include "Core.h"

#include "FormDouble.h"
#include "Pars.h"
#include "ParsStochasticTrajectory.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

using trajectory::ParsStochastic;
using trajectory::ParsEvolved;
using std::string;

namespace pythonext {

PARS_GETTER_SETTER(unsigned long, ParsStochastic, seed)
PARS_GETTER_SETTER(bool, ParsStochastic, noise)
PARS_GETTER_SETTER(size_t, ParsStochastic, nTraj)


void export_15_ParsStochasticTrajectory()
{
  scope namespaceScope = trajectoryNameSpace;
  class_<ParsStochastic, bases<ParsEvolved> >
    (
      "ParsStochastic",
      "Wrapper of :core:`trajectory::ParsStochastic`",
      init<parameters::ParameterTable&, optional<const std::string&> >()
        [with_custodian_and_ward<1,2>()]
    )
    .PARS_PROPERTY(seed)
    .PARS_PROPERTY(noise)
    .PARS_PROPERTY(nTraj)
  ;
}

} // pythonext