// -*- C++ -*-

#include "PythonExtension.h"

#include <StateVector.tcc>

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/stringize.hpp>

using namespace boost::python;

using quantumdata::StateVector;


namespace pythonext{

void export_StateVector()
{

  // TODO: add some useful constructors, e.g. allow to initialize a StateVector from a numpy array

#define DECL_DIRECT_PRODUCTS(z,r2,data) .def(self * other< StateVector<BOOST_PP_ADD(r2,1)> >())
#define BOOST_PP_LOCAL_MACRO(n) class_<StateVector<n> >(BOOST_PP_STRINGIZE(BOOST_PP_CAT(StateVector, n)), no_init) \
  .def(self + self) \
  .def(self - self) \
  .def(self * dcomp()) \
  .def(dcomp() * self) \
  .def(self / dcomp()) \
  .def("norm", &StateVector<n>::norm) \
  .def("renorm", &StateVector<n>::renorm) \
  BOOST_PP_REPEAT(BOOST_PP_SUB(PYTHON_MAX_RANK,n), DECL_DIRECT_PRODUCTS, n) ;
#define BOOST_PP_LOCAL_LIMITS (1, PYTHON_MAX_RANK)
#include BOOST_PP_LOCAL_ITERATE()

}

} // pythonext

