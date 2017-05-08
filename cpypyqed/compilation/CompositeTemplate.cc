// Copyright Raimar Sandner 2012–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "Composite.tcc"
#include "Act.h"
#include "QuantumSystem.h"

#include <boost/python.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/shared_ptr.hpp>

#define CPPQED_CORE_GIT "@CPPQED_CORE_GIT@"

using namespace boost::python;
using composite::result_of::Make;

namespace pythonext {{

namespace {{
typedef Make<{actlist} >::type ThisCompositeConstPtr;
typedef boost::remove_const<ThisCompositeConstPtr::element_type>::type ThisComposite;
static const int RANK=composite::MaxRank<ThisComposite::Acts>::type::value+1;

const ThisCompositeConstPtr (*ptrToMake)({const_act_arguments}) = &composite::make;
}}

BOOST_PYTHON_MODULE({modulename})
{{
  object obj=class_<ThisComposite,bases<structure::QuantumSystem<RANK> > >("{classname}",no_init);
  obj.attr("core_git") = std::string(CPPQED_CORE_GIT);
  def("compositeMake",  ptrToMake,
    {custodian_and_wards}
  );
  register_ptr_to_python<ThisCompositeConstPtr>();
}}

}} // pythonext