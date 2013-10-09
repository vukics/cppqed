// -*- C++ -*-
#ifndef   UTILS_EVOLVEDGSL_TCC_INCLUDED
#define   UTILS_EVOLVEDGSL_TCC_INCLUDED

#include "EvolvedGSL.h"

// #include "ArrayTraits.h"
// If we include this here, without having the system see also the necessary overload-specializations, then we get an error at linking time because the general definition is never found.
// Instead, we rely on inclusion-order dependence: the necessary overload-declarations must be included before EvolvedGSL.tcc
#include "Evolved.tcc"

#include <boost/make_shared.hpp>


namespace evolved {

  
namespace details  {

// Two indirections needed because cannot declare member functions for Impl here

ImplPtr createImpl(void*, size_t, int(double,const double*,double*,void*), double, double, const double*, SteppingFunction);

void apply(ImplPtr, double*, double, double*, double*);

size_t extractFailedSteps(ImplPtr);

} // details


template<typename A>
class MakerGSL<A>::GSL : public Evolved<A> 
{
public:
  typedef Evolved<A> Base;

  typedef typename Base::Derivs Derivs;

  using Base::getA; using Base::getTime; using Base::getDtTry;

  GSL(A& a, Derivs derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs, SteppingFunction sf, double nextDtTryCorretionFactor)
    : Base(a,derivs,dtInit,epsRel,epsAbs),
      pImpl_(details::createImpl(this,cpputils::size(getA()),auxFunction,epsRel,epsAbs,cpputils::data(scaleAbs),sf)),
      sf_(sf),
      nextDtTryCorretionFactor_(nextDtTryCorretionFactor)
  {
    if (!cpputils::isStorageContiguous(a)) throw (NonContiguousStorageException());
  }

  std::ostream& displayParameters_v(std::ostream& os) const {return os<<"# EvolvedGSL implementation, stepping function: "<<sf_<<std::endl;}

private:

  static int auxFunction(double t, const double* y, double* dydt, void* aux)
  {
    using namespace cpputils;
    
    GSL* e=static_cast<GSL*>(aux);

    const A    yInterfaceA(create(y   ,e->getA()));
          A dydtInterfaceA(create(dydt,e->getA()));

    e->getDerivs()(t,yInterfaceA,dydtInterfaceA);

    return details::onSuccess;
  }

  void step_v(double deltaT)
  {
    double
      time=getTime(),
      dtTry=getDtTry(),
      nextDtTry=( fabs(deltaT)<fabs(dtTry/nextDtTryCorretionFactor_) ? dtTry : 0. );

    apply(pImpl_,&time,time+deltaT,&dtTry,cpputils::data(getA()));

    Base::update(time, nextDtTry ? nextDtTry/nextDtTryCorretionFactor_ : dtTry );

  }

  size_t nFailedSteps_v() const {return extractFailedSteps(pImpl_);}

  const details::ImplPtr pImpl_;

  const SteppingFunction sf_;
  const double nextDtTryCorretionFactor_;

};




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
  return boost::make_shared<GSL, A& >(a,derivs,dtInit,epsRel,epsAbs,scaleAbs,sf_,nextDtTryCorretionFactor_);
}



} // evolved

#endif // UTILS_EVOLVEDGSL_TCC_INCLUDED
