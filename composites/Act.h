// -*- C++ -*-
#ifndef   _ACT_HELPER_INCLUDED
#define   _ACT_HELPER_INCLUDED

#include "SubSystem.h"

#include "TMP_Tools.h"

#include <boost/mpl/size.hpp>


#define BASE_class1 tmptools::Vector<BOOST_PP_ENUM_PARAMS(TMPTOOLS_MAX_VECTOR_SIZE,V)>
#define BASE_class2 structure::SubSystemsInteraction<mpl::size<BASE_class1>::value>

template<BOOST_PP_ENUM_BINARY_PARAMS(TMPTOOLS_MAX_VECTOR_SIZE,int V,=tmptools::vectorDefaultArgument BOOST_PP_INTERCEPT)>
class Act
  : public BASE_class1,
    public BASE_class2
{
public:
  typedef          BASE_class1              Vector     ;
  typedef typename BASE_class2::Interaction Interaction;

  explicit Act(const Interaction& ia) : BASE_class2(&ia) {}

};

#undef  BASE_class2
#undef  BASE_class1


#endif // _ACT_HELPER_INCLUDED
