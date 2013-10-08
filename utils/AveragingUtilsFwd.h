// -*- C++ -*-
#ifndef   ELEMENTS_UTILS_AVERAGINGUTILSFWD_H_INCLUDED
#define   ELEMENTS_UTILS_AVERAGINGUTILSFWD_H_INCLUDED


template<int>
class ReducedDensityOperator;

namespace averagingUtils {
  
template<int, bool IS_TIME_DEPENDENT=false>
class Collecting;

template<int, int, bool IS_TIME_DEPENDENT>
class Transferring;

} // averagingUtils


#endif // ELEMENTS_UTILS_AVERAGINGUTILSFWD_H_INCLUDED
