// Copyright Raimar Sandner 2012–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-

#include "Act.h"
#include "Interaction.h"

#include <string>

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

#define CPPQED_CORE_GIT "@CPPQED_CORE_GIT@"

using namespace boost::python;

namespace pythonext {{

BOOST_PYTHON_MODULE({modulename})
{{
  object obj=class_<{act} >
    (
      "{classname}",
     init<const {act}::InteractionPtr::element_type&>()
        [with_custodian_and_ward<1,2>()]
    )
  ;
  obj.attr("core_git") = std::string(CPPQED_CORE_GIT);

}}

}} // pythonext