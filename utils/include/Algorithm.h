// -*- C++ -*-
#ifndef   _CPPUTILS_ALGORITHM_H
#define   _CPPUTILS_ALGORITHM_H

namespace cpputils {


template<typename In, typename In2, typename BinOp>
BinOp 
for_each(In, In, In2, BinOp);


template<typename In, typename T, typename UnOp>
T
accumulate(In, In, T, UnOp);


template<typename In, typename T, typename UnOp, typename BinOp>
T
accumulate(In, In, T, UnOp, BinOp);


template<typename SeqOfSeqs, typename Out>
const Out&
concatenate(const SeqOfSeqs&, Out&);


} // cpputils


#include "impl/Algorithm.tcc"

#endif // _CPPUTILS_ALGORITHM_H
