
#include "utils.h"

using namespace boost::python;

namespace pythonext {

const PyArrayObject * numeric_np(const numeric::array &arr)
{
  if(!PyArray_Check(arr.ptr())){
    PyErr_SetString(PyExc_ValueError, "expected a PyArrayObject");
    throw_error_already_set();
  }
  return reinterpret_cast<const PyArrayObject *>(arr.ptr());
}

} // pythonext