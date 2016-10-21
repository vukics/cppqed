// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMDATA_NEGPT_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMDATA_NEGPT_TCC_INCLUDED

#include "core_config.h"

#include "NegPT.h"

#ifndef DO_NOT_USE_FLENS

#include "Blitz2FLENS.tcc"
#include "BlitzArraySliceIterator.tcc"


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

  blitzplusplus::basi::Transposer<2*RANK,ExtendV>::transpose(rhoShallowPT);

  DensityOperatorLow rhoDeepPT(rhoShallowPT.shape()); rhoDeepPT=rhoShallowPT;

  return sum(blitzplusplus::selectNegative(real(blitz2flens::ev(rhoDeepPT))));

}


} // quantumdata

#endif // DO_NOT_USE_FLENS

#endif // CPPQEDCORE_QUANTUMDATA_NEGPT_TCC_INCLUDED
