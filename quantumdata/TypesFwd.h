// -*- C++ -*-
#ifndef   _STRUCTURE_TYPES_FWD_H
#define   _STRUCTURE_TYPES_FWD_H


namespace structure {

struct StaticTag;
// We declare this here because this is the only file which is included everywhere in structure.

} // structure


namespace quantumdata {

namespace details {

struct Empty;

} // details

template<int, typename=details::Empty>
struct Types;

} // quantumdata

#endif // _STRUCTURE_TYPES_FWD_H
