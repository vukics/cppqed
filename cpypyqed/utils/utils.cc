
#include "utils.h"

#include <boost/lexical_cast.hpp>
#include <string>

using namespace boost::python;

namespace pythonext {

const PyArrayObject * numeric_np(const numeric::array &arr, size_t rank)
{
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  const PyArrayObject *result = reinterpret_cast<const PyArrayObject *>(arr.ptr());
  if(rank && rank!=PyArray_NDIM(const_cast<PyArrayObject *>(result))) {
    PyErr_SetString(PyExc_RuntimeError, (std::string("Expected an array with rank ")+boost::lexical_cast<std::string>(rank)).c_str());
    throw_error_already_set();
  }
  return result;
}

} // pythonext