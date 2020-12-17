// Copyright Raimar Sandner 2012–2019.
// Copyright András Vukics 2019–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#ifndef CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
#define CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED

#include "BlitzTiny.h"
#include "BlitzArray.h"
#include "ComplexExtensions.h"

#include <boost/python/numpy.hpp>

namespace pythonext {

template<typename T, int RANK>
auto arrayToNumpy(const blitz::Array<T,RANK>& a)
{
  using namespace boost::python::numpy;
  std::vector<long> shape(a.shape().begin(),a.shape().end()), strides(a.stride().begin(),a.stride().end());
  return from_data(a.copy().dataFirst(),dtype::get_builtin<T>(),shape,strides,boost::python::object());
}

template<typename T, int RANK>
auto numpyToArray(const boost::python::numpy::ndarray& array)
{
  using std::string; using namespace boost::python;
  ExtTiny<RANK> shape, strides;
  for (size_t i=0; i<RANK; ++i) {
    shape(i)=array.shape(i);
    strides(i)=array.strides(i);
  }
  return blitz::Array<T,RANK>(reinterpret_cast<T*>(array.get_data()),shape,strides,blitz::duplicateData);
}

} // pythonext

#endif // CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
