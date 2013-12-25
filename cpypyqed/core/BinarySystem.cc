// -*- C++ -*-

#include "PythonExtension.h"
#include "Namespaces.h"

#include "BinarySystem.h"
#include "Interaction.h"

#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/stringize.hpp>

using namespace boost::python;

using structure::QuantumSystem;

namespace pythonext {

void export_QuantumSystem_BinarySystem()
{
  {
    scope namespaceScope = binaryNameSpace;
    class_<binary::Base, bases<QuantumSystem<2> >, boost::noncopyable >("Base", no_init);
    def("make", binary::doMake,
        with_custodian_and_ward_postcall<0, 1>()
    );
    binaryNameSpace.staticmethod("make");
  }
#define BOOLS (true)(false)
#define BINARYSYSTEM_INSTANTIATIONS(z,b) \
  class_<BinarySystem<BOOST_PP_SEQ_ENUM(b)>, bases<binary::Base>, boost::noncopyable> \
      (BOOST_PP_STRINGIZE( \
        BOOST_PP_CAT( \
          BOOST_PP_CAT( \
            BOOST_PP_CAT( \
              BOOST_PP_CAT( \
                BOOST_PP_CAT(BinarySystem_,BOOST_PP_SEQ_ELEM(0,b)) \
                ,_) \
              ,BOOST_PP_SEQ_ELEM(1,b)) \
            ,_) \
          ,BOOST_PP_SEQ_ELEM(2,b)) \
        ),no_init);
BOOST_PP_SEQ_FOR_EACH_PRODUCT(BINARYSYSTEM_INSTANTIATIONS, (BOOLS)(BOOLS)(BOOLS))

  register_ptr_to_python<boost::shared_ptr<binary::Base const> >();
}

} //pythonext