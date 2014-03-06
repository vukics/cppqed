// -*- C++ -*-

#ifndef CPYPYQED_INCLUDE_UTILS_H_INCLUDED
#define CPYPYQED_INCLUDE_UTILS_H_INCLUDED

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "PythonExtension.h"

#include <numpy/ndarrayobject.h>

namespace pythonext {

const PyArrayObject * numeric_np(const boost::python::numeric::array &arr);

} //pythonext

#endif // CPYPYQED_INCLUDE_UTILS_H_INCLUDED