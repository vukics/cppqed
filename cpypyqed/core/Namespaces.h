// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

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
class py_binary{};
class py_evolution{};

extern boost::python::class_<py_trajectory>        trajectoryNameSpace;
extern boost::python::class_<py_evolved>           evolvedNameSpace;
extern boost::python::class_<py_structure>         structureNameSpace;
extern boost::python::class_<py_quantumdata>       quantumdataNameSpace;
extern boost::python::class_<py_parameters>        parametersNameSpace;
extern boost::python::class_<py_quantumtrajectory> quantumtrajectoryNameSpace;
extern boost::python::class_<py_binary>            binaryNameSpace;
extern boost::python::class_<py_evolution>         evolutionNameSpace;

} // pythonext

#endif // CPYPYQED_CORE_NAMESPACES_H_INCLUDED