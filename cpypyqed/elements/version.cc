// Copyright Raimar Sandner 2012–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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