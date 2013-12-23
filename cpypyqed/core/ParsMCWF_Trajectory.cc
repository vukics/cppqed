// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Pars.h"
#include "ParsMCWF_Trajectory.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

using quantumtrajectory::ParsMCWF;
using trajectory::ParsStochastic;
using std::string;

namespace pythonext {

PARS_GETTER_SETTER(double,   ParsMCWF, dpLimit)
PARS_GETTER_SETTER(double,   ParsMCWF, overshootTolerance)
PARS_GETTER_SETTER(int,      ParsMCWF, logLevel)

void export_20_ParsMCWF_Trajectory()
{
  scope namespaceScope = quantumtrajectoryNameSpace;
  class_<ParsMCWF, bases<ParsStochastic> >("ParsMCWF",init<parameters::ParameterTable&, optional<const std::string&> >())
    .PARS_PROPERTY(dpLimit)
    .PARS_PROPERTY(overshootTolerance)
    .PARS_PROPERTY(logLevel)
  ;
}

} // pythonext