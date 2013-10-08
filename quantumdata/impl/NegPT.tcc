// -*- C++ -*-
#ifndef   QUANTUMDATA_IMPL_NEGPT_TCC_INCLUDED
#define   QUANTUMDATA_IMPL_NEGPT_TCC_INCLUDED

#include "NegPT.h"

#ifndef DO_NOT_USE_FLENS

#include "impl/Blitz2FLENS.tcc"
#include "impl/BlitzArraySliceIterator.tcc"

#include <boost/mpl/transform.hpp>

namespace quantumdata {


namespace details {


inline double selectNegative(double d) {return d<0 ? d : 0;}

BZ_DECLARE_FUNCTION_RET(selectNegative,double)

template<int, typename>
struct ExtendV;


} // details


template<int RANK, typename V>
double negPT(const DensityOperator<RANK>& rho, V)
{
  using namespace blitz2flens;

  typedef typename DensityOperator<RANK>::DensityOperatorLow DensityOperatorLow;

  typedef typename details::ExtendV<RANK,V>::type ExtendedV;

  typedef GeMatrixMF<dcomp,RowMajor>::type GeMatrix;

  DensityOperatorLow rhoShallowPT(rho());

  blitzplusplus::basi::Transposer<2*RANK,ExtendedV>::transpose(rhoShallowPT);

  DensityOperatorLow rhoDeepPT(rhoShallowPT.shape()); rhoDeepPT=rhoShallowPT;

  CArray<1> eigenValues(rho.getTotalDimension());

  {
    GeMatrix a(matrix<RowMajor>(rhoDeepPT));
    typename DenseVectorMF<dcomp>::type v(blitz2flens::vector(eigenValues));

    ev(false,false,a,v,a,a);
  }

  return sum(details::selectNegative(real(eigenValues)));

}


namespace details {


namespace namehider {

using namespace boost::mpl;
namespace mpl=boost::mpl;
using namespace tmptools;

template<int RANK, typename V>
struct Algorithm 
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
	 >
{};

} // namehider


template<int RANK, typename V>
struct ExtendV : namehider::Algorithm<RANK,V>
{};


} // details


} // quantumdata

#endif // DO_NOT_USE_FLENS

#endif // QUANTUMDATA_IMPL_NEGPT_TCC_INCLUDED