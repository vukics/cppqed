// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Core.h"

#include "Pars.h"

#include <boost/python/docstring_options.hpp>

using namespace boost::python;

namespace pythonext {

void update(parameters::Table &p, list args, str s=str("--")){
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
  class_<parameters::Table,boost::noncopyable>("ParameterTable", "Wrapper of :core:`parameters::Table`")
    .def("printList", &parameters::Table::printList)
    ;
  docOpts.disable_cpp_signatures();
  def("update",update,
R"doc(Wrapper of :core:`parameters::update`. Note that the signature differs a little bit
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