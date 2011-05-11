// -*- C++ -*-
#ifndef   STRUCTURE_ELEMENT_LIOUVILLEAN_FWD_INCLUDED
#define   STRUCTURE_ELEMENT_LIOUVILLEAN_FWD_INCLUDED

namespace structure {

template<int RANK, int NOJ=1, bool IS_TD=false> // NOJ stands for the number of jumps
// Note that even an elementary system can have several possible
// jumps. Eg direction of recoil for atoms. This should be known at
// compile time.
class ElementLiouvillean;


} // structure

#endif // STRUCTURE_ELEMENT_LIOUVILLEAN_FWD_INCLUDED
