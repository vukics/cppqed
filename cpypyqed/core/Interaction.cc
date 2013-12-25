// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Interaction.h"

using namespace boost::python;

using structure::Interaction;

namespace pythonext {

void export_Interaction() {
  scope namespaceScope = structureNameSpace;
  class_<Interaction<2>, boost::noncopyable >("Interaction2", no_init);
  class_<Interaction<3>, boost::noncopyable >("Interaction3", no_init);
}


} // pythonext