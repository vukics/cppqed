// -*- C++ -*-

#include "Composite.tcc"
#include "Act.h"
#include "QuantumSystem.h"

#include <boost/python.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost::python;
using composite::result_of::Make;

namespace pythonext {{

namespace {{
typedef Make<Act<0,1>,Act<0,2> >::type ThisCompositeConstPtr;
typedef boost::remove_const<ThisCompositeConstPtr::element_type>::type ThisComposite;
static const int RANK=composite::MaxRank<ThisComposite::Acts>::type::value+1;

const ThisCompositeConstPtr (*ptrToMake)({const_act_arguments}) = &composite::make;
}}

BOOST_PYTHON_MODULE({modulename})
{{
  class_<ThisComposite,bases<structure::QuantumSystem<RANK> > >("{classname}",no_init);
  def("compositeMake",  ptrToMake,
    {custodian_and_wards}
  );
  register_ptr_to_python<ThisCompositeConstPtr>();
}}

}} // pythonext