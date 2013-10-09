// -*- C++ -*-
#ifndef   ELEMENTS_COMPOSITES_ACT_H_INCLUDED
#define   ELEMENTS_COMPOSITES_ACT_H_INCLUDED

#include "SubSystem.h"

#include "SmartPtr.h"
#include "TMP_Tools.h"

#include <boost/mpl/size.hpp>


#define BASE_class1 tmptools::Vector<V...>
#define BASE_class2 composite::SubSystemsInteraction<mpl::size<BASE_class1>::value>

template<int... V>
class Act
  : public BASE_class1,
    public BASE_class2
{
public:
  typedef BASE_class1 Vector;

  template<typename IA>
  explicit Act(const IA& ia) : BASE_class2(cpputils::sharedPointerize(ia)) {}

};

#undef  BASE_class2
#undef  BASE_class1


#endif // ELEMENTS_COMPOSITES_ACT_H_INCLUDED
