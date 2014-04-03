// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "JaynesCummings_.h"
#include "Mode.tcc"
#include "Pars.h"
#include "Qbit_.h"

#include "ParsPropertyMacros.h"

using namespace boost::python;

using jaynescummings::Pars;
using jaynescummings::Base;

namespace pythonext{

PARS_GETTER_SETTER(dcomp, Pars, g)

namespace {
const jaynescummings::Ptr jc_make(const qbit::Ptr& q, const mode::Ptr& m)
{
  return jaynescummings::make(q,m,1);
}
}

void export_JaynesCummings()
{

  class_<Base<false>, boost::noncopyable> ("JaynesCummingsBase", no_init);
  {
    scope namespaceScope = jaynescummingsNameSpace;
    class_<Pars>
      (
        "Pars",
      init<parameters::ParameterTable&, optional<const std::string&> >()
      [with_custodian_and_ward<1, 2>()]
      )
      .PARS_PROPERTY(g)
    ;
    register_ptr_to_python< jaynescummings::Ptr >();
    implicitly_convertible<boost::shared_ptr<Base<false>>,jaynescummings::Ptr>();

    def("make", jaynescummings::make<qbit::Ptr, mode::Ptr>,
        with_custodian_and_ward_postcall<0,1,with_custodian_and_ward_postcall<0,2,with_custodian_and_ward_postcall<0,3>>>());
    jaynescummingsNameSpace.staticmethod("make");

    implicitly_convertible<Pars, dcomp>();
    implicitly_convertible<boost::shared_ptr<Base<false>>, boost::shared_ptr<const structure::Interaction<2>>>();
  }
}

} // pythonext