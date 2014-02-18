/// \briefFile{Defines the helper class composite::_}
#ifndef   COMPOSITES_ACT_H_INCLUDED
#define   COMPOSITES_ACT_H_INCLUDED

#include "SubSystem.h"

#include "SmartPtr.h"
#include "TMP_Tools.h"

#include <boost/mpl/size.hpp>


/// Auxiliary tools to Composite
namespace composite {


#define BASE_class SubSystemsInteraction<mpl::size<tmptools::Vector<V...> >::value>

template<int... V>
class _
  : public tmptools::Vector<V...>,
    public BASE_class
{
public:
  typedef tmptools::Vector<V...> Vector;

  template<typename IA>
  explicit _(const IA& ia) : BASE_class(cpputils::sharedPointerize(ia)) {}

};

#undef  BASE_class


}


/// Template alias for backward compatibility
template<int... V>
using Act = composite::_<V...>;


#endif // COMPOSITES_ACT_H_INCLUDED
