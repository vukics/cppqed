// -*- C++ -*-
#ifndef   STRUCTURE_AVERAGED_FWD_INCLUDED
#define   STRUCTURE_AVERAGED_FWD_INCLUDED

namespace structure {

class AveragedCommon;

#ifndef   NDEBUG
struct AveragesFishyException;
#endif // NDEBUG

template<int, bool IS_TD=true> class Averaged;

} // structure


#endif // STRUCTURE_AVERAGED_FWD_INCLUDED
