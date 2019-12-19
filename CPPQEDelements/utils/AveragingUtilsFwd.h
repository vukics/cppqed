// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_UTILS_AVERAGINGUTILSFWD_H_INCLUDED
#define   CPPQEDELEMENTS_UTILS_AVERAGINGUTILSFWD_H_INCLUDED


template<int>
class ReducedDensityOperator;

namespace averagingUtils {
  
template<int, bool IS_TIME_DEPENDENT=false>
class Collecting;

template<int, int, bool IS_TIME_DEPENDENT>
class Transferring;

} // averagingUtils


#endif // CPPQEDELEMENTS_UTILS_AVERAGINGUTILSFWD_H_INCLUDED
