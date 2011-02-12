// -*- C++ -*-
#ifndef   STRUCTURE_ELEMENT_LIOUVILLEAN_IMPL_INCLUDED
#define   STRUCTURE_ELEMENT_LIOUVILLEAN_IMPL_INCLUDED

#include "Range.h"

#include<boost/bind.hpp>



namespace structure {


template<int RANK, int NOJ>
const LiouvilleanCommon::Probabilities ElementLiouvillean<RANK,NOJ>::probabilities(const LazyDensityOperator& matrix) const
{
  Probabilities probas(NOJ);
  // Note that this cannot be anything like static because of the by reference semantics of blitz::Array

  boost::transform(jumpProbas_,probas.begin(),
		   bind(&JumpProbabilityStrategy::operator(),_1,boost::cref(matrix)));
  return probas;
}



template<int RANK>
const LiouvilleanCommon::Probabilities ElementLiouvillean<RANK,1>::probabilities(const LazyDensityOperator& M) const
{
  Probabilities probas(1); probas(0)=probability(M); return probas;
}


} // structure


#endif // STRUCTURE_ELEMENT_LIOUVILLEAN_IMPL_INCLUDED
