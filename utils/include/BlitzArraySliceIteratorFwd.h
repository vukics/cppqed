// -*- C++ -*-
#ifndef   _BLITZ_ARRAY_SMART_ITERATOR_FWD_H
#define   _BLITZ_ARRAY_SMART_ITERATOR_FWD_H


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



#endif // _BLITZ_ARRAY_SMART_ITERATOR_FWD_H
