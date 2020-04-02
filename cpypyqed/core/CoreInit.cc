// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "PythonExtension.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL core_ARRAY_API
#include <numpy/arrayobject.h>

namespace pythonext {

void export_00_Init()
{
  import_array1();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
}

} // pythonext
