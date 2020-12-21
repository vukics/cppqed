// Copyright Raimar Sandner 2012–2019.
// Copyright András Vukics 2019–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#ifndef CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
#define CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED

#include "BlitzTiny.h"

#include <boost/python/numpy.hpp>

namespace cpputils {

template<typename T, int RANK>
auto arrayToNumpy(const blitz::Array<T,RANK>& a)
{
  using namespace boost::python::numpy;
  Py_intptr_t shape[RANK]; // Plain array is the best here: nothing has the temptation on the C++ side to delete it (apparently, it happens on the Python side …)
  std::copy(a.shape().begin(),a.shape().end(),&shape[0]);
  auto res{empty(RANK,shape,dtype::get_builtin<T>())}; // The cleanest solution: owns its data by default. But, it necessitates the below raw-pointer solution.
  std::copy(a.begin(),a.end(),reinterpret_cast<T*>(res.get_data()));
  return res;
}

template<typename T, int RANK>
auto numpyToArray(const boost::python::numpy::ndarray& array)
{
  using namespace boost::python;
  ExtTiny<RANK> shape, strides;
  for (size_t i=0; i<RANK; ++i) {
    shape(i)=array.shape(i);
    strides(i)=array.strides(i);
  }
  return blitz::Array<T,RANK>(reinterpret_cast<T*>(array.get_data()),shape,strides,blitz::duplicateData);
}

} // cpputils

#endif // CPYPYQED_INCLUDE_BLITZ2NUMPY_TCC_INCLUDED
