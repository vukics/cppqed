// -*- C++ -*-
#ifndef   UTILS_RANGE_H_INCLUDED
#define   UTILS_RANGE_H_INCLUDED

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



#endif // UTILS_RANGE_H_INCLUDED
