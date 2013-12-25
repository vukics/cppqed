// -*- C++ -*-

#ifndef CPYPYQED_ELEMENTS_NAMESPACES_H_INCLUDED
#define CPYPYQED_ELEMENTS_NAMESPACES_H_INCLUDED

#include "PythonExtension.h"

namespace pythonext {

class py_particle{};
class py_qbit{};
class py_mode{};
class py_particlecavity{};
class py_jaynescummings{};

extern boost::python::class_<py_particle>       particleNameSpace;
extern boost::python::class_<py_qbit>           qbitNameSpace;
extern boost::python::class_<py_mode>           modeNameSpace;
extern boost::python::class_<py_particlecavity> particlecavityNameSpace;
extern boost::python::class_<py_jaynescummings> jaynescummingsNameSpace;

} // pythonext

#endif // CPYPYQED_ELEMENTS_NAMESPACES_H_INCLUDED