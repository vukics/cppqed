// -*- C++ -*-
// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_DETAILS_BINARYINTERACTIONGENERATOR_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_DETAILS_BINARYINTERACTIONGENERATOR_H_INCLUDED

#include "SmartPtr.h"


struct EmptyAveragingBaseForBinaryInteractions {};

#endif // ELEMENTS_INTERACTIONS_DETAILS_BINARYINTERACTIONGENERATOR_H_INCLUDED


template<typename A=EmptyAveragingBaseForBinaryInteractions>
class BIG_CLASS_NAME : public BIG_NAMESPACE_NAME::Base, public A
{
public:

  template<typename F1, typename F2>
  BIG_CLASS_NAME(const F1& f1, const F2& f2 BIG_ADDITIONAL_PARAMETERS)
    : BIG_NAMESPACE_NAME::Base(cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f2) BIG_ADDITIONAL_PARAMETERS_PASS),
      A()
  {}

};



#undef BIG_ADDITIONAL_PARAMETERS_PASS
#undef BIG_ADDITIONAL_PARAMETERS
#undef BIG_NAMESPACE_NAME
#undef BIG_CLASS_NAME
