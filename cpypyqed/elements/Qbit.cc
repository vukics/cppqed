// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Qbit_.h"
#include "ParsQbit.h"
#include "QM_PictureFwd.h"

using namespace boost::python;

using quantumdata::StateVector;

namespace pythonext {

void export_Qbit()
{
  class_<QbitBase, bases<structure::QuantumSystem<1>>, boost::noncopyable>("QbitBase",no_init);
  {
    scope namespaceScope = qbitNameSpace;
    def("make", qbit::make, with_custodian_and_ward_postcall<0,1>());
    qbitNameSpace.staticmethod("make");
    register_ptr_to_python< qbit::Ptr >();
    implicitly_convertible<boost::shared_ptr<QbitBase>,qbit::Ptr>();

    def("state0", qbit::state0); qbitNameSpace.staticmethod("state0");
    def("state1", qbit::state1); qbitNameSpace.staticmethod("state1");
    def("init", (const StateVector<1> (*)(const dcomp&)) &qbit::init);
    def("init", (const StateVector<1> (*)(const qbit::Pars&)) &qbit::init);
    qbitNameSpace.staticmethod("init");
  }
}

  
} // pythonext