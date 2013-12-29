// -*- C++ -*-
#include "PythonExtension.h"

#include "elements_version.h"

#include <string>

using namespace boost::python;

namespace pythonext {

void export_version()
{
  scope currentScope;
  currentScope.attr("elements_git") = std::string(g_CPPQEDelements_GIT_SHA1);
}

} // pythonext