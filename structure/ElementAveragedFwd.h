// -*- C++ -*-
#ifndef   STRUCTURE_ELEMENT_AVERAGED_FWD_INCLUDED
#define   STRUCTURE_ELEMENT_AVERAGED_FWD_INCLUDED

namespace structure {

class ElementAveragedCommon;

template<int RANK, bool IS_TD=false>
class ElementAveraged;

template<int RANK, bool IS_TD=false>
class ClonableElementAveraged;


namespace averaged {

class DiagonalDO;

} // averaged

} // structure

#endif // STRUCTURE_ELEMENT_AVERAGED_FWD_INCLUDED
