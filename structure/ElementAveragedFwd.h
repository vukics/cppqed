// -*- C++ -*-
#ifndef   STRUCTURE_ELEMENTAVERAGEDFWD_H_INCLUDED
#define   STRUCTURE_ELEMENTAVERAGEDFWD_H_INCLUDED

namespace structure {


template<int RANK, bool IS_TD=false>
class ElementAveraged;

template<int RANK, bool IS_TD=false>
class ClonableElementAveraged;


namespace averaged {

class ReducedDensityOperator;

template<int, bool IS_TD=false>
class Collecting;

} // averaged

} // structure

#endif // STRUCTURE_ELEMENTAVERAGEDFWD_H_INCLUDED
