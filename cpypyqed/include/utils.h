// -*- C++ -*-

#ifndef CPYPYQED_INCLUDE_UTILS_H_INCLUDED
#define CPYPYQED_INCLUDE_UTILS_H_INCLUDED

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "PythonExtension.h"

#include <numpy/arrayobject.h>

namespace pythonext {

inline const PyArrayObject * numeric_np(const boost::python::numeric::array &arr, size_t rank=0)
{
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    boost::python::throw_error_already_set();
  }
  PyArrayObject *result = reinterpret_cast<PyArrayObject *>(arr.ptr()); // cannot be const because of NPY API weirdness
  if(rank && rank!=PyArray_NDIM(result)) {
    PyErr_SetString(PyExc_RuntimeError, (std::string("Expected an array with rank ")+boost::lexical_cast<std::string>(rank)).c_str());
    boost::python::throw_error_already_set();
  }
  return result;
}

} //pythonext

#endif // CPYPYQED_INCLUDE_UTILS_H_INCLUDED