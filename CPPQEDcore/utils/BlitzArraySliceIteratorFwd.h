// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BLITZARRAYSLICEITERATORFWD_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZARRAYSLICEITERATORFWD_H_INCLUDED


namespace blitzplusplus {

namespace basi {


template<int RANK, typename V, bool IS_CONST>
class Iterator;

template<int RANK, typename V>
class Transposer;

template<int RANK, typename V>
class Indexer;


} // basi


template<int RANK, typename V>
class SlicesData;


namespace basi_fast {


template<int RANK, typename V, bool IS_CONST>
class Iterator;


} // basi_fast


} // blitzplusplus



#endif // CPPQEDCORE_UTILS_BLITZARRAYSLICEITERATORFWD_H_INCLUDED
