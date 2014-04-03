// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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