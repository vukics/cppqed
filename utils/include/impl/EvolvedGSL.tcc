// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_EVOLVEDGSL_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_EVOLVEDGSL_TCC_INCLUDED

#include "EvolvedGSL.h"

#include "impl/Evolved.tcc"

#include <boost/make_shared.hpp>

namespace evolved {


template<typename A>
const typename Maker<A>::Ptr MakerGSL<A>::operator()(
							  A& a,
							  Derivs derivs,
							  double dtInit,
							  double epsRel,
							  double epsAbs,
							  const A& scaleAbs
							  ) const
{
  return boost::make_shared<details::EvolvedGSL<A>, A& >(a,derivs,dtInit,epsRel,epsAbs,scaleAbs,sf_,nextDtTryCorretionFactor_);
}



namespace details {

template<typename A>
int auxFunction(double t, const double* y, double* dydt, void* aux)
{
  EvolvedGSL<A>* e=static_cast<EvolvedGSL<A>*>(aux);

  const A    yInterfaceA(cpputils::ArrayMemoryTraits<A>::create(y   ,e->getA()));
        A dydtInterfaceA(cpputils::ArrayMemoryTraits<A>::create(dydt,e->getA()));

  e->getDerivs()(t,yInterfaceA,dydtInterfaceA);

  return onSuccess;
}



template<typename A>
EvolvedGSL<A>::EvolvedGSL(
			  A& a,
			  Derivs derivs,
			  double dtInit,
			  double epsRel,
			  double epsAbs,
			  const A& scaleAbs,
			  SteppingFunction sf,
			  double nextDtTryCorretionFactor
			  )
  : Base(a,derivs,dtInit,epsRel,epsAbs),
    pImpl_(createImpl(this,Traits::size(getA()),auxFunction<A>,epsRel,epsAbs,Traits::data(scaleAbs),sf)),
    sf_(sf),
    nextDtTryCorretionFactor_(nextDtTryCorretionFactor)
{
  if (!Traits::isStorageContiguous(a)) throw (NonContiguousStorageException());
}



template<typename A> void EvolvedGSL<A>::doStep(double deltaT)
{
  double
    time=getTime(),
    dtTry=getDtTry(),
    nextDtTry=( fabs(deltaT)<fabs(dtTry/nextDtTryCorretionFactor_) ? dtTry : 0. );

  apply(pImpl_,&time,time+deltaT,&dtTry,Traits::data(getA()));

  Base::update(time, nextDtTry ? nextDtTry/nextDtTryCorretionFactor_ : dtTry );

}


} // details


} // evolved

#endif // UTILS_INCLUDE_IMPL_EVOLVEDGSL_TCC_INCLUDED
