// -*- C++ -*-

#include "DispatchFreeType.h"


template<typename A=EmptyAveragingBaseForBinaryInteractions>
class BIG_CLASS_NAME : public BIG_NAMESPACE_NAME::Base, public A
{
public:

  template<typename F1, typename F2>
  BIG_CLASS_NAME(const F1& f1, const F2& f2 BIG_ADDITIONAL_PARAMETERS)
    : BIG_NAMESPACE_NAME::Base(dispatchFreeType(f1),dispatchFreeType(f2) BIG_ADDITIONAL_PARAMETERS_PASS),
      A()
  {}

};



#undef BIG_ADDITIONAL_PARAMETERS_PASS
#undef BIG_ADDITIONAL_PARAMETERS
#undef BIG_NAMESPACE_NAME
#undef BIG_CLASS_NAME
