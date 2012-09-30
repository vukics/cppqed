// -*- C++ -*-
#ifndef   STRUCTURE_AVERAGEDFWD_H_INCLUDED
#define   STRUCTURE_AVERAGEDFWD_H_INCLUDED

namespace structure {

class AveragedCommon;

#ifndef   NDEBUG
struct AveragesNumberMismatchException;
#endif // NDEBUG

template<int, bool IS_TD=true> class Averaged;

} // structure


#endif // STRUCTURE_AVERAGEDFWD_H_INCLUDED
