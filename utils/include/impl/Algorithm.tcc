// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_ALORITHM_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_ALORITHM_TCC_INCLUDED

#include "Algorithm.h"

#include<boost/bind.hpp>

#include<functional>
#include<numeric>


namespace cpputils {


template<typename In, typename In2, typename BinOp>
BinOp 
for_each(In first, In last, In2 first2, BinOp op)
{
  while (first!=last) op(*first++,*first2++);
  return op;
} // NEEDS_WORK does BOOST.zip_iterator obsolate this?


template<typename In, typename T, typename UnOp>
T
accumulate(In first, In last, T init, UnOp op)
{
  return std::accumulate(first,last,init,bind(std::plus<T>(),_1,bind(op,_2)));
}

template<typename In, typename T, typename UnOp, typename BinOp>
T
accumulate(In first, In last, T init, UnOp op, BinOp binOp)
{
  return std::accumulate(first,last,init,bind(binOp,_1,bind(op,_2)));
}



namespace details {

template<typename Iter, typename T>
inline
const Iter
concatenateHelper(Iter iter, const T& t)
{
  return std::copy(t.begin(),t.end(),iter);
}

} // details


template<typename SeqOfSeqs_In, typename Out_Iterator>
const Out_Iterator
concatenate(SeqOfSeqs_In begin, SeqOfSeqs_In end, Out_Iterator out)
{
  return std::accumulate(begin,end,out,details::concatenateHelper<Out_Iterator,typename SeqOfSeqs_In::value_type>); 
  // concatenation is accumulation (!)
}


template<typename SeqOfSeqs, typename Out>
const Out
concatenateGrow(const SeqOfSeqs& sOs, Out empty)
{
  concatenate(sOs.begin(),sOs.end(),std::back_inserter(empty));
  return empty;
}


template<typename SeqOfSeqs, typename Out>
const Out&
concatenate(const SeqOfSeqs& sOs, Out& out)
{
  concatenate(sOs.begin(),sOs.end(),out.begin());
  return out;
}



} // cpputils


#endif // UTILS_INCLUDE_IMPL_ALORITHM_TCC_INCLUDED
