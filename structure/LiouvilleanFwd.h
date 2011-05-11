// -*- C++ -*-
#ifndef   STRUCTURE_LIOUVILLEAN_FWD_INCLUDED
#define   STRUCTURE_LIOUVILLEAN_FWD_INCLUDED


namespace structure {


struct LiouvilleanCommon;

#ifndef   NDEBUG
struct LiouvilleanFishyException;
#endif // NDEBUG

template<int, bool IS_TD=true>
class Liouvillean;


} // structure


#endif // STRUCTURE_LIOUVILLEAN_FWD_INCLUDED
