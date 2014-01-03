#include "PythonExtension.h"
#include "Namespaces.h"

using namespace boost::python;

namespace pythonext {

class_<py_trajectory>         trajectoryNameSpace("trajectory");
class_<py_evolved>            evolvedNameSpace("evolved");
class_<py_structure>          structureNameSpace("structure");
class_<py_quantumdata>        quantumdataNameSpace("quantumdata");
class_<py_parameters>         parametersNameSpace("parameters");
class_<py_quantumtrajectory>  quantumtrajectoryNameSpace("quantumtrajectory");
class_<py_binary>             binaryNameSpace("binary");

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
}

} // pythonext