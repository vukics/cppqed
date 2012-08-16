// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_MULTIINDEXITERATOR_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_MULTIINDEXITERATOR_TCC_INCLUDED

#include "MultiIndexIterator.h"

namespace cpputils {

namespace details {


template<int RANK>
void doIt(const TTD_IDXTINY(RANK)&, 
	  const TTD_IDXTINY(RANK)&,
	  TTD_IDXTINY(RANK)& idx,
	  boost::mpl::int_<0>)
{
  idx(0)++;
  // This will of course eventually put the iterator into an illegal
  // state when idx(0)>ubound(0), but this is how every (unchecked)
  // iterator works.
}



template<int RANK, int N>
void doIt(const TTD_IDXTINY(RANK)& lbound, 
	  const TTD_IDXTINY(RANK)& ubound,
	  TTD_IDXTINY(RANK)& idx,
	  boost::mpl::int_<N>)
{
  if (idx(N)==ubound(N)) {
    idx(N)=lbound(N);
    doIt(lbound,ubound,idx,boost::mpl::int_<N-1>());
  }
  else idx(N)++;
}



} // details

} // cpputils


#endif // UTILS_INCLUDE_IMPL_MULTIINDEXITERATOR_TCC_INCLUDED
