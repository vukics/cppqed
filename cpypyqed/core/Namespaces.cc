// Copyright Raimar Sandner 2012–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Core.h"
#include "Namespaces.h"

using namespace boost::python;

namespace pythonext {

class_<py_trajectory>         trajectoryNameSpace("trajectory", "The :core:`trajectory` namespace");
class_<py_evolved>            evolvedNameSpace("evolved", "The :core:`evolved` namespace");
class_<py_structure>          structureNameSpace("structure", "The :core:`structure` namespace");
class_<py_quantumdata>        quantumdataNameSpace("quantumdata", "The :core:`quantumdata` namespace");
class_<py_parameters>         parametersNameSpace("parameters", "The :core:`parameters` namespace");
class_<py_quantumtrajectory>  quantumtrajectoryNameSpace("quantumtrajectory", "The :core:`quantumtrajectory` namespace");
class_<py_binary>             binaryNameSpace("binary", "The :core:`binary` namespace");
class_<py_evolution>          evolutionNameSpace("evolution", "The :core:`evolution` namespace");

void export_0Namespaces()
{
  scope current;
  current.attr("trajectory")        = trajectoryNameSpace;
  current.attr("evolved")           = evolvedNameSpace;
  current.attr("structure")         = structureNameSpace;
  current.attr("quantumdata")       = quantumdataNameSpace;
  current.attr("parameters")        = parametersNameSpace;
  current.attr("quantumtrajectory") = quantumtrajectoryNameSpace;
  current.attr("binary")            = binaryNameSpace;
  current.attr("evolution")         = evolutionNameSpace;
}

} // pythonext