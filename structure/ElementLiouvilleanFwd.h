// -*- C++ -*-
#ifndef   STRUCTURE_ELEMENT_LIOUVILLEAN_FWD_INCLUDED
#define   STRUCTURE_ELEMENT_LIOUVILLEAN_FWD_INCLUDED

namespace structure {

template<int RANK, int NOJ=1> // NOJ stands for the number of jumps
// Note that even an elementary system can have several possible
// jumps. Eg direction of recoil for atoms. This should be known at
// compile time.
class ElementLiouvillean;

template<int RANK>
class ElementLiouvillean<RANK,1>;

} // structure

#endif // STRUCTURE_ELEMENT_LIOUVILLEAN_FWD_INCLUDED
