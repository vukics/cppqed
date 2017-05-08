// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the helper class composite::_}
#ifndef   CPPQEDCORE_COMPOSITES_ACT_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_ACT_H_INCLUDED

#include "SubSystem.h"

#include "SmartPtr.h"
#include "TMP_Tools.h"

#include <boost/mpl/size.hpp>


/// Auxiliary tools to Composite
namespace composite {


#define BASE_class SubSystemsInteraction<mpl::size<tmptools::Vector<V...> >::value>

/// Helper class to Composite
/**
 * It combines a set of \link retainedindexpositionsdefined retained index positions\endlink (via its tmptools::Vector base) with a certain structure::Interaction element
 * (via the SubSystemsInteraction base) of the corresponding arity (equalling the number of retained index positions) – this means that the given structure::Interaction
 * element will act on the given retained index positions.
 *
 * Composite expects a \refBoost{Boost.Fusion,fusion} sequence as its template argument `VA`.
 *
 * \tparam V has the same role as in tmptools::Vector
 */
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


#endif // CPPQEDCORE_COMPOSITES_ACT_H_INCLUDED
