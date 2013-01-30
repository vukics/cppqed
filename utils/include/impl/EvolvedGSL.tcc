// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_EVOLVEDGSL_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_EVOLVEDGSL_TCC_INCLUDED

#include "EvolvedGSL.h"

// #include "ArrayTraits.h"
// If we include this here, without having the system see also the necessary overload-specializations, then we get an error at linking time because the general definition is never found.
// Instead, we rely on inclusion-order dependence: the necessary overload-declarations must be included before EvolvedGSL.tcc
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
  return boost::make_shared<details::GSL<A>, A& >(a,derivs,dtInit,epsRel,epsAbs,scaleAbs,sf_,nextDtTryCorretionFactor_);
}



namespace details {

using namespace cpputils;
  
template<typename A>
int auxFunction(double t, const double* y, double* dydt, void* aux)
{
  GSL<A>* e=static_cast<GSL<A>*>(aux);

  const A    yInterfaceA(create(y   ,e->getA()));
        A dydtInterfaceA(create(dydt,e->getA()));

  e->getDerivs()(t,yInterfaceA,dydtInterfaceA);

  return onSuccess;
}



template<typename A>
GSL<A>::GSL(
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
    pImpl_(createImpl(this,cpputils::size(getA()),auxFunction<A>,epsRel,epsAbs,cpputils::data(scaleAbs),sf)),
    sf_(sf),
    nextDtTryCorretionFactor_(nextDtTryCorretionFactor)
{
  if (!cpputils::isStorageContiguous(a)) throw (NonContiguousStorageException());
}



template<typename A> void GSL<A>::doStep(double deltaT)
{
  double
    time=getTime(),
    dtTry=getDtTry(),
    nextDtTry=( fabs(deltaT)<fabs(dtTry/nextDtTryCorretionFactor_) ? dtTry : 0. );

  apply(pImpl_,&time,time+deltaT,&dtTry,cpputils::data(getA()));

  Base::update(time, nextDtTry ? nextDtTry/nextDtTryCorretionFactor_ : dtTry );

}


} // details


} // evolved

#endif // UTILS_INCLUDE_IMPL_EVOLVEDGSL_TCC_INCLUDED
