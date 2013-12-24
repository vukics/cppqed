// -*- C++ -*-

#ifndef CPYPYQED_CORE_NAMESPACES_H_INCLUDED
#define CPYPYQED_CORE_NAMESPACES_H_INCLUDED

#include "PythonExtension.h"

namespace pythonext {

class py_trajectory{};
class py_evolved{};
class py_structure{};
class py_quantumdata{};
class py_parameters{};
class py_quantumtrajectory{};

extern boost::python::class_<py_trajectory>        trajectoryNameSpace;
extern boost::python::class_<py_evolved>           evolvedNameSpace;
extern boost::python::class_<py_structure>         structureNameSpace;
extern boost::python::class_<py_quantumdata>       quantumdataNameSpace;
extern boost::python::class_<py_parameters>        parametersNameSpace;
extern boost::python::class_<py_quantumtrajectory> quantumtrajectoryNameSpace;

} // pythonext

#endif // CPYPYQED_CORE_NAMESPACES_H_INCLUDED