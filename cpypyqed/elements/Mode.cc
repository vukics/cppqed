// Copyright Raimar Sandner 2012–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)

#include "PythonExtension.h"
#include "Namespaces.h"

#include "Mode.tcc"
#include "ParsMode.h"
#include "QM_PictureFwd.h"
#include "StateVector.tcc"

using namespace boost::python;

namespace pythonext {
//BOOST_PYTHON_FUNCTION_OVERLOADS(make_overloads,mode::make,2,3)

namespace {
const mode::Ptr (*make1)(const mode::Pars&, QM_Picture) = &mode::make;
const mode::Ptr (*make2)(const mode::ParsPumped&, QM_Picture) = &mode::make;
const mode::Ptr (*make3)(const mode::ParsLossy&, QM_Picture) = &mode::make;
const mode::Ptr (*make4)(const mode::ParsPumpedLossy&, QM_Picture) = &mode::make;
}

void export_Mode()
{

  class_<ModeBase, bases<structure::QuantumSystem<1>>, boost::noncopyable>("ModeBase",no_init);

  {
    scope namespaceScope = modeNameSpace;
    def("make", make1, with_custodian_and_ward_postcall<0, 1>());
    def("make", make2, with_custodian_and_ward_postcall<0, 1>());
    def("make", make3, with_custodian_and_ward_postcall<0, 1>());
    def("make", make4, with_custodian_and_ward_postcall<0, 1>());
    modeNameSpace.staticmethod("make");
    register_ptr_to_python< mode::Ptr >();
    implicitly_convertible<boost::shared_ptr<ModeBase>, mode::Ptr>();

    def("coherent", mode::coherent); modeNameSpace.staticmethod("coherent");
    def("fock", mode::fock);         modeNameSpace.staticmethod("fock");
    def("init", mode::init);         modeNameSpace.staticmethod("init");
  }
}

} // pythonext