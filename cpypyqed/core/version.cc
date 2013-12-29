// -*- C++ -*-
#include "PythonExtension.h"

#include "core_version.h"

#include <string>

using namespace boost::python;

namespace pythonext {

void export_version()
{
  scope currentScope;
  currentScope.attr("core_git") = std::string(g_CPPQEDcore_GIT_SHA1);
}

} // pythonext