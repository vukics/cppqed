// Copyright Raimar Sandner 2012–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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