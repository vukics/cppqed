// Copyright Raimar Sandner 2012–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "PythonExtension.h"
#include "Namespaces.h"

using namespace boost::python;

namespace pythonext {

class_<py_particle>       particleNameSpace("particle");
class_<py_qbit>           qbitNameSpace("qbit");
class_<py_mode>           modeNameSpace("mode");
class_<py_particlecavity> particlecavityNameSpace("particlecavity");
class_<py_jaynescummings> jaynescummingsNameSpace("jaynescummings");

void export_0Namespaces()
{
  scope current;
  current.attr("particle")       = particleNameSpace;
  current.attr("qbit")           = qbitNameSpace;
  current.attr("mode")           = modeNameSpace;
  current.attr("particlecavity") = particlecavityNameSpace;
  current.attr("jaynescummings") = jaynescummingsNameSpace;
}

} // pythonext