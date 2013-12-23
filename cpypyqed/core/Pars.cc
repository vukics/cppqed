// -*- C++ -*-
#include "PythonExtension.h"
#include "Namespaces.h"

#include "Pars.h"

using namespace boost::python;

namespace pythonext {

void update(parameters::ParameterTable &p, list args, str s=str("--")){
  int argc = len(args);
  char **argv = new char*[argc];
  for (int c=0; c<argc; c++) {
    argv[c] = extract<char *>(args[c]);
  }
  parameters::update(p,argc,argv,extract<std::string>(s));
  delete[] argv;
}


void export_05_Pars()
{
  scope namespaceScope = parametersNameSpace;
  class_<parameters::ParameterTable,boost::noncopyable>("ParameterTable")
    .def("printList", &parameters::ParameterTable::printList)
    ;
  def("update",update);
}

} // pythonext