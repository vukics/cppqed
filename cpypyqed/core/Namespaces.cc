#include "PythonExtension.h"
#include "Namespaces.h"

using namespace boost::python;

namespace pythonext {

object trajectoryNameSpace  = class_<py_trajectory>    ("trajectory");
object evolvedNameSpace     = class_<py_evolved>       ("evolved");
object structureNameSpace   = class_<py_structure>     ("structure");
object quantumdataNameSpace = class_<py_quantumdata>   ("quantumdata");
object parametersNameSpace  = class_<py_parameters>    ("paramters");

void export_0Namespaces()
{
  scope current;
  current.attr("trajectory")  = trajectoryNameSpace;
  current.attr("evolved")     = evolvedNameSpace;
  current.attr("structure")   = structureNameSpace;
  current.attr("quantumdata") = quantumdataNameSpace;
  current.attr("parameters")  = parametersNameSpace;
}

} // pythonext