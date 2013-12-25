// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "QuantumSystem.h"

#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/stringize.hpp>

using namespace boost::python;

using structure::QuantumSystem;


namespace pythonext {

void export_QuantumSystem(){

  scope namespaceScope = structureNameSpace;
#define QUANTUMSYSTEM_INSTANTIATIONS(z,r,data) \
  class_<QuantumSystem<r>, boost::noncopyable >(BOOST_PP_STRINGIZE(BOOST_PP_CAT(QuantumSystem,r)), no_init) \
    .def("highestFrequency",&QuantumSystem<r>::highestFrequency) \
  ;
BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_ADD(PYTHON_HALF_RANK,1), QUANTUMSYSTEM_INSTANTIATIONS, data)

}

} // pythonext