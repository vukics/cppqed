// -*- C++ -*-

#include "Act.h"
#include "Interaction.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost::python;

namespace pythonext {{

BOOST_PYTHON_MODULE({modulename})
{{
  class_<{act} >
    (
      "{classname}",
     init<const {act}::InteractionPtr::element_type&>()
        [with_custodian_and_ward<1,2>()]
    )
  ;
}}

}} // pythonext