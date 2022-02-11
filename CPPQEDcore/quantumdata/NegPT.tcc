// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMDATA_NEGPT_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_NEGPT_TCC_INCLUDED

#include "NegPT.h"

#ifndef DO_NOT_USE_FLENS

#include "Blitz2FLENS.h"
#include "SliceIterator.tcc"


namespace quantumdata {


template<int RANK, typename V>
double negPT(const DensityOperator<RANK>& rho, V)
{
  using namespace boost::mpl;
  namespace mpl=boost::mpl;
  using namespace tmptools;

  struct ExtendV
    : fold<Range<RANK,RANK>,
          typename fold<Ordinals<RANK>,
                        vector_c<int>,
                        push_back<mpl::_1,
                                  if_<numerical_contains<V,mpl::_2>,
                                      plus<mpl::_2,int_<RANK> >,
                                      mpl::_2
                                      >
                                  >
                        >::type,
          push_back<mpl::_1,
                    if_<numerical_contains<V,minus<mpl::_2,int_<RANK> > >,
                        minus<mpl::_2,int_<RANK> >,
                        mpl::_2
                        >
                    >
          >::type
  {};

  typedef typename DensityOperator<RANK>::DensityOperatorLow DensityOperatorLow;

  DensityOperatorLow rhoShallowPT(rho.getArray());

  cppqedutils::sliceiterator::Transposer<CArray,2*RANK,ExtendV>::_(rhoShallowPT);

  DensityOperatorLow rhoDeepPT(rhoShallowPT.shape()); rhoDeepPT=rhoShallowPT;

  return sum(blitzplusplus::selectNegative(real(blitz2flens::ev(rhoDeepPT))));

}


} // quantumdata

#endif // DO_NOT_USE_FLENS

#endif // CPPQEDCORE_QUANTUMDATA_NEGPT_TCC_INCLUDED
