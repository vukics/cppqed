// -*- C++ -*-
#ifndef   STRUCTURE_ELEMENTLIOUVILLEANFWD_H_INCLUDED
#define   STRUCTURE_ELEMENTLIOUVILLEANFWD_H_INCLUDED

namespace structure {

template<int RANK, int NOJ=1, bool IS_TIME_DEPENDENT=false> // NOJ stands for the number of jumps
// Note that even an elementary system can have several possible
// jumps. Eg direction of recoil for atoms. This should be known at
// compile time.
class ElementLiouvillean;

template<int RANK, int NOJ, bool IS_TIME_DEPENDENT=false>
class ElementLiouvilleanStrategies;

} // structure

#endif // STRUCTURE_ELEMENTLIOUVILLEANFWD_H_INCLUDED