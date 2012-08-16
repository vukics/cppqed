// -*- C++ -*-
#ifndef   QUANTUMDATA_TYPESFWD_H_INCLUDED
#define   QUANTUMDATA_TYPESFWD_H_INCLUDED


namespace structure {

struct StaticTag;
// We declare this here because this is the only file which is included everywhere in structure.

} // structure


namespace quantumdata {

namespace details {

struct Empty;

} // details

template<int RANK, typename B=details::Empty>
struct Types;

} // quantumdata

#endif // QUANTUMDATA_TYPESFWD_H_INCLUDED
