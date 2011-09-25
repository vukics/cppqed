// -*- C++ -*-
#ifndef   _ACT_HELPER_INCLUDED
#define   _ACT_HELPER_INCLUDED

#include "SubSystem.h"

#include "TMP_Tools.h"

#include<boost/preprocessor/repetition.hpp>

#include<boost/mpl/size.hpp>

#define VECTOR_MAX_SIZE 20

#define BASE_class1 tmptools::Vector<BOOST_PP_ENUM_PARAMS(VECTOR_MAX_SIZE,V)>
#define BASE_class2 structure::SubSystemsInteraction<mpl::size<BASE_class1>::value>

template<BOOST_PP_ENUM_BINARY_PARAMS(VECTOR_MAX_SIZE,int V,=TMPTOOLS_VECTOR_DEFAULT_ARG BOOST_PP_INTERCEPT)>
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

#undef  VECTOR_MAX_SIZE

/*
#include<boost/preprocessor/iteration/iterate.hpp>
#include<boost/preprocessor/repetition.hpp>
#include<boost/preprocessor/arithmetic/sub.hpp>

#define BOOST_PP_ITERATION_LIMITS (2,10)
#define BOOST_PP_FILENAME_1 "details/ActImplementationSpecializations.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS
*/

#endif // _ACT_HELPER_INCLUDED
