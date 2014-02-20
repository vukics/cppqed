// -*- C++ -*-

#ifndef CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
#define CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "PythonExtension.h"

#include "BlitzTiny.h"
#include "ComplexExtensions.h"

#include <boost/type_traits/is_same.hpp>
#include <numpy/ndarrayobject.h>

namespace pythonext {

template<typename A, int RANK>
boost::python::object arrayToNumpy(const A &a)
{
  using namespace boost::python;

  ExtTiny<RANK> dims;
  npy_intp npy_dims[RANK];
  dims = a.extent();
  std::copy(dims.begin(),dims.end(), npy_dims);
  PyObject * pyObj;
  if(boost::is_same<dcomp,typename A::T_numtype>::value) {
    pyObj = PyArray_SimpleNewFromData(RANK, npy_dims, NPY_CDOUBLE, const_cast<A&>(a).dataFirst());
  }
  else if(boost::is_same<double,typename A::T_numtype>::value) {
    pyObj = PyArray_SimpleNewFromData(RANK, npy_dims, NPY_DOUBLE, const_cast<typename A::T_numtype *>(a.dataFirst()));
  }
  else {
    PyErr_SetString(PyExc_NotImplementedError, "Conversion of a blitz array with this data type not implemented.");
    throw_error_already_set();
  }
  handle<> h( pyObj );
  numeric::array arr( h );
  return arr.copy();
}

} // pythonext

#endif // CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED