// -*- C++ -*-

#include "PythonExtension.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_ARRAY_UNIQUE_SYMBOL core_ARRAY_API
#include <numpy/arrayobject.h>

namespace pythonext {

void export_00_Init()
{
  import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
}

} // pythonext