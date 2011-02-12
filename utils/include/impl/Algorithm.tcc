// -*- C++ -*-
#ifndef   _CPPUTILS_ALGORITHM_IMPL_H
#define   _CPPUTILS_ALGORITHM_IMPL_H

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
} // NEEDS_WORK does BOOST.zip_iterator eliminate the need for this?


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


template<typename SeqOfSeqs, typename Out>
const Out&
concatenate(const SeqOfSeqs& sOs, Out& out)
{
  std::accumulate(sOs.begin(),sOs.end(),out.begin(),details::concatenateHelper<typename Out::iterator,typename SeqOfSeqs::value_type>); 
  // concatenation is accumulation (!)
  return out;
}



} // cpputils


#endif // _CPPUTILS_ALGORITHM_IMPL_H
