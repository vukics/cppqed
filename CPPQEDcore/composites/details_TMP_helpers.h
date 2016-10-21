// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_COMPOSITES_DETAILS_TMP_HELPERS_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_DETAILS_TMP_HELPERS_H_INCLUDED

#include "TMP_Tools.h"

#include <boost/fusion/include/fold.hpp>
#include <boost/mpl/less.hpp>

namespace mpl=boost::mpl;

namespace composite {

using namespace mpl;

template<typename VEC, typename Pred=less<mpl::_1,mpl::_2> >
struct MaxMF : deref<typename boost::mpl::max_element<VEC,Pred>::type>::type {};

template<typename VEC1, typename VEC2>
struct SeqLess : less<typename MaxMF<VEC1>::type,typename MaxMF<VEC2>::type> {};

namespace namehider {

using namespace tmptools;

template<typename VA, typename ICW>
struct InnerCheck : fold<VA,false_,or_<numerical_contains<mpl::_2,ICW>,mpl::_1> > {};

template<int RANK, typename VA>
struct Algorithm
  : fold<Ordinals<RANK>,
         true_,
         and_<mpl::_1,
              InnerCheck<VA,mpl::_2>
              >
         > 
{};


} // namehider


template<int RANK, typename VA>
struct CheckMeta : namehider::Algorithm<RANK,VA>
{};


} // composite


#endif // CPPQEDCORE_COMPOSITES_DETAILS_TMP_HELPERS_H_INCLUDED
