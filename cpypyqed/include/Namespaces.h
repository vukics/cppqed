// -*- C++ -*-

#ifndef CPYPYQED_INCLUDE_NAMESPACES_H_INCLUDED
#define CPYPYQED_INCLUDE_NAMESPACES_H_INCLUDED

#include "PythonExtension.h"

namespace pythonext {

class py_trajectory{};
class py_evolved{};
class py_structure{};
class py_quantumdata{};
class py_parameters{};

extern boost::python::object trajectoryNameSpace;
extern boost::python::object evolvedNameSpace;
extern boost::python::object structureNameSpace;
extern boost::python::object quantumdataNameSpace;
extern boost::python::object parametersNameSpace;

} // pythonext

#endif // CPYPYQED_INCLUDE_NAMESPACES_H_INCLUDED