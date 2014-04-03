// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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
