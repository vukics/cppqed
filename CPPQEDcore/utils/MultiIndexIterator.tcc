// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_MULTIINDEXITERATOR_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_MULTIINDEXITERATOR_TCC_INCLUDED

#include "MultiIndexIterator.h"

namespace cpputils {

template<int RANK>
void MultiIndexIterator<RANK>::doIt(boost::mpl::int_<0>)
{
  idx_(0)++;
  // This will of course eventually put the iterator into an illegal state when idx(0)>ubound(0), but this is how every (unchecked) iterator works.
}



template<int RANK> template<int N>
void MultiIndexIterator<RANK>::doIt(boost::mpl::int_<N>)
{
  if (idx_(N)==ubound_(N)) {
    idx_(N)=lbound_(N);
    doIt(boost::mpl::int_<N-1>());
  }
  else idx_(N)++;
}

} // cpputils


#endif // CPPQEDCORE_UTILS_MULTIINDEXITERATOR_TCC_INCLUDED
