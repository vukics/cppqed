// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "Core.h"

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
  class_<QuantumSystem<r>, boost::noncopyable >(BOOST_PP_STRINGIZE(BOOST_PP_CAT(QuantumSystem,r)),\
                                                "Instantiation of :core:`structure::QuantumSystem` with RANK="#r,\
                                                no_init) \
    .def("highestFrequency",&QuantumSystem<r>::highestFrequency) \
  ; \
  register_ptr_to_python<QuantumSystem<r>::Ptr>(); \
  implicitly_convertible<boost::shared_ptr<QuantumSystem<r>>, QuantumSystem<r>::Ptr>();
BOOST_PP_REPEAT_FROM_TO(1, BOOST_PP_ADD(PYTHON_HALF_RANK,1), QUANTUMSYSTEM_INSTANTIATIONS, data)

}

} // pythonext