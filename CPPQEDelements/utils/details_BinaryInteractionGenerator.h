// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{A “preprocessor factory” to binary interactions}
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_UTILS_DETAILS_BINARYINTERACTIONGENERATOR_H_INCLUDED
#define   CPPQEDELEMENTS_UTILS_DETAILS_BINARYINTERACTIONGENERATOR_H_INCLUDED

#include "SmartPtr.h"

struct EmptyAveragingBaseForInteractions {};

#endif // CPPQEDELEMENTS_UTILS_DETAILS_BINARYINTERACTIONGENERATOR_H_INCLUDED


#ifndef   BIG_ADDITIONAL_TEMPLATE_PARAMETERS
#define   BIG_ADDITIONAL_TEMPLATE_PARAMETERS
#endif // BIG_ADDITIONAL_TEMPLATE_PARAMETERS

#ifndef   BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS
#define   BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS
#endif // BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS

#ifndef   BIG_ADDITIONAL_PARAMETERS
#define   BIG_ADDITIONAL_PARAMETERS
#endif // BIG_ADDITIONAL_PARAMETERS

#ifndef   BIG_ADDITIONAL_PARAMETERS_PASS
#define   BIG_ADDITIONAL_PARAMETERS_PASS
#endif // BIG_ADDITIONAL_PARAMETERS_PASS


#define BASE_class BIG_NAMESPACE_NAME::Base BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS


template<BIG_ADDITIONAL_TEMPLATE_PARAMETERS typename AveragingType=EmptyAveragingBaseForInteractions>
class BIG_CLASS_NAME : public BASE_class, public AveragingType
{
public:

  template<typename F1, typename F2, typename... AveragingConstructorParameters>
  BIG_CLASS_NAME(const F1& f1, const F2& f2 BIG_ADDITIONAL_PARAMETERS , AveragingConstructorParameters&&... a)
    : BASE_class(cpputils::sharedPointerize(f1),cpputils::sharedPointerize(f2) BIG_ADDITIONAL_PARAMETERS_PASS),
      AveragingType(std::forward<AveragingConstructorParameters>(a)...)
  {}

};


#undef BASE_class

#undef BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS
#undef BIG_ADDITIONAL_TEMPLATE_PARAMETERS
#undef BIG_ADDITIONAL_PARAMETERS_PASS
#undef BIG_ADDITIONAL_PARAMETERS
#undef BIG_NAMESPACE_NAME
#undef BIG_CLASS_NAME
