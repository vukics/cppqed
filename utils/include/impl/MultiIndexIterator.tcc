// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_MULTIINDEXITERATOR_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_MULTIINDEXITERATOR_TCC_INCLUDED

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


#endif // UTILS_INCLUDE_IMPL_MULTIINDEXITERATOR_TCC_INCLUDED
