// -*- C++ -*-
#include "PythonExtension.h"
#include "Namespaces.h"

#include "Pars.h"

#include <boost/python/docstring_options.hpp>

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
  docstring_options docOpts;
  scope namespaceScope = parametersNameSpace;
  class_<parameters::ParameterTable,boost::noncopyable>("ParameterTable", "Wrapper of :core:`parameters::ParameterTable`")
    .def("printList", &parameters::ParameterTable::printList)
    ;
  docOpts.disable_cpp_signatures();
  def("update",update,
R"doc(Wrapper of :core:`parameters::ParameterTable::update`. Note that the signature differs a little bit
compared to the C++ version, it is not necessary to pass in the argument list length.

:param p: The :class:`ParameterTable`
:type p: :class:`ParameterTable`
:param list argv: The argument list
:param str prefix: The prefix (default "``--``"))doc",
       (arg("p"),"argv",arg("prefix")="--")
     );
  parametersNameSpace.staticmethod("update");
}

} // pythonext