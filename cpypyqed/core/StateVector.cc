// -*- C++ -*-

#include "Core.h"
#include "blitz2numpy.tcc"

#include <StateVector.tcc>

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <boost/python/import.hpp>

using namespace boost::python;

using quantumdata::StateVector;


namespace pythonext{

template<int RANK>
struct SV_to_python_numpy
{
  static PyObject* convert(const StateVector<RANK> &s)
  {
    object sv = import(BOOST_PP_STRINGIZE(BOOST_PP_CAT(cpypyqed, DEBUG_SUFFIX)) ".pycppqed.statevector");
    object StateVector = sv.attr("StateVector");
    return boost::python::incref(StateVector(arrayToNumpy<CArray<RANK>,RANK>(s.getArray())).ptr());
  }
};


void export_StateVector()
{

  // TODO: add some useful constructors, e.g. allow to initialize a StateVector from a numpy array

  {
    scope namespaceScope = quantumdataNameSpace;

#define DECL_DIRECT_PRODUCTS(z,r2,data) .def(self * other< StateVector<BOOST_PP_ADD(r2,1)> >())
#define BOOST_PP_LOCAL_MACRO(n) to_python_converter<StateVector<n>, SV_to_python_numpy<n>>();
#define BOOST_PP_LOCAL_LIMITS (1, PYTHON_HALF_RANK)
#include BOOST_PP_LOCAL_ITERATE()

  }
}

} // pythonext

