// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#ifndef CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
#define CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "PythonExtension.h"

#include "BlitzTiny.h"
#include "BlitzArray.h"
#include "ComplexExtensions.h"
#include "utils.h"

#include <blitz/array.h>
#include <boost/type_traits/is_same.hpp>
#include <boost/lexical_cast.hpp>
#include <numpy/ndarrayobject.h>

namespace pythonext {

namespace {

template<typename A>
NPY_TYPES npy_type(const A&)
{
  using namespace boost::python;
  if(cpputils::TypeID<A>::value == "DArray") return NPY_DOUBLE;
  if(cpputils::TypeID<A>::value == "CArray") return NPY_CDOUBLE;
  PyErr_SetString(PyExc_NotImplementedError, "Conversion of a blitz array with this data type not implemented.");
  throw_error_already_set();
  return NPY_NOTYPE;
}

} // anonymous namespace


template<int RANK>
ExtTiny<RANK> numeric_shape(const boost::python::numeric::array &a)
{
  const PyArrayObject *np_a=numeric_np(a, RANK);
  npy_intp *dims=PyArray_DIMS(const_cast<PyArrayObject *>(np_a));
  ExtTiny<RANK> shape;
  std::copy(dims,dims+RANK,shape.begin());
  return shape;
}


template<typename A, int RANK>
boost::python::object arrayToNumpy(const A &a)
{
  using namespace boost::python;

  ExtTiny<RANK> dims;
  npy_intp npy_dims[RANK];
  dims = a.extent();
  std::copy(dims.begin(),dims.end(), npy_dims);
  PyObject * pyObj = PyArray_SimpleNewFromData(RANK, npy_dims, npy_type(a), const_cast<A&>(a).dataFirst());;
  handle<> h( pyObj );
  numeric::array arr( h );
  return arr.copy();
}

template<typename dtype, int RANK>
blitz::Array<dtype,RANK> numpyToArray(const boost::python::numeric::array &array)
{
  using std::string;
  using namespace boost::python;
  NPY_TYPES type = npy_type(blitz::Array<dtype,RANK>(0));
  ExtTiny<RANK> shape = numeric_shape<RANK>(array);
  PyArrayObject *casted = (PyArrayObject*)PyArray_FromAny(array.ptr(),PyArray_DescrFromType(type), 0, PYTHON_MAX_RANK, 0, NULL);
  blitz::Array<dtype,RANK> blitz_a = blitz::Array<dtype,RANK>(static_cast<dtype *>(PyArray_DATA(casted)), shape, blitz::duplicateData);
  return blitz_a;
}

} // pythonext

#endif // CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
