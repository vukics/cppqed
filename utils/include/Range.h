// -*- C++ -*-
#ifndef   _RANGE_UTILS_INCLUDED
#define   _RANGE_UTILS_INCLUDED

#include "Algorithm.h"

#include <boost/range.hpp>
#include <boost/range/algorithm.hpp>


namespace cpputils {


template<typename Range, typename In2, class Function>
Function
for_each(const Range& r, In2 i2, Function f)
{
  return for_each(boost::begin(r),boost::end(r),i2,f);
}

template<typename Range, typename In2, class Function>
Function
for_each(      Range& r, In2 i2, Function f)
{
  return for_each(boost::begin(r),boost::end(r),i2,f);
}


} // cpputils



#endif // _RANGE_UTILS_INCLUDED
