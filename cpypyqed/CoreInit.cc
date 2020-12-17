// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "PythonExtension.h"

namespace pythonext {

void export_00_Init()
{
  Py_Initialize();
  boost::python::numpy::initialize();
}

} // pythonext
