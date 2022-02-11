// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_COMPOSITES_DETAILS_TMP_HELPERS_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_DETAILS_TMP_HELPERS_H_INCLUDED

#include "TMP_Tools.h"

#include <boost/fusion/include/fold.hpp>
#include <boost/mpl/less.hpp>

namespace mpl=boost::mpl;

namespace composite {

using namespace mpl;

template<typename VEC, typename Pred=less<mpl::_1,mpl::_2> >
using MaxMF_t = typename deref<typename boost::mpl::max_element<VEC,Pred>::type>::type;

template<typename VEC1, typename VEC2>
struct SeqLess : less<MaxMF_t<VEC1>,MaxMF_t<VEC2>> {};

namespace namehider {

template<typename VA, typename ICW>
struct InnerCheck : fold<VA,false_,or_<tmptools::numerical_contains<mpl::_2,ICW>,mpl::_1> > {};

} // namehider


template<int RANK, typename VA>
constexpr bool CheckMeta_v=
  fold<tmptools::Ordinals<RANK>,
       true_,
       and_<mpl::_1,
            namehider::InnerCheck<VA,mpl::_2>
           >
      >::type::value;


} // composite


#endif // CPPQEDCORE_COMPOSITES_DETAILS_TMP_HELPERS_H_INCLUDED
