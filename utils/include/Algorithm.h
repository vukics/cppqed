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


template<typename SeqOfSeqs_In, typename Out_Iterator>
const Out_Iterator
concatenate(SeqOfSeqs_In begin, SeqOfSeqs_In end, Out_Iterator);


template<typename SeqOfSeqs, typename Out>
const Out
concatenateGrow(const SeqOfSeqs&, Out empty);
// Filling an empty container with concatenated values


template<typename SeqOfSeqs, typename Out>
const Out&
concatenate(const SeqOfSeqs&, Out&);
// Filling a container of the necessary size with concatenated values


} // cpputils


#endif // _CPPUTILS_ALGORITHM_H
