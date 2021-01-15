// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Implementation of MakerGSL::_}
// Note that some ArrayTraits file declaring the used traits templates must be included *before* inclusion of this file.
// If the traits are also instantiated after inclusion of this file, then the ArrayTraits file included beforehand must also contain the traits *definitions*.
#ifndef   CPPQEDCORE_UTILS_EVOLVEDGSL_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_EVOLVEDGSL_TCC_INCLUDED

#include "EvolvedGSL.h"


namespace evolved {

  
namespace details { // We declare indirections here because at this point we cannot yet declare member functions for Impl
                    // (as Impl at this point is forward-declared only, and will be defined only in EvolvedGSL.cc according to the rationale that the framework’s header files must not contain GSL headers).

ImplPtr createImpl(void*, size_t, int(double,const double*,double*,void*), double, double, const double*, SteppingFunction);

void apply(ImplPtr, double*, double, double*, double*);

size_t extractFailedSteps(ImplPtr);

} // details


/// The class that actually implements the Evolved interface
/**
 * It bridges the low-level \GSL and higher levels that work with the array type `A`. For this, it has to connect in a bi-directional way the array type `A` with C-arrays used by \GSL,
 * which is performed by the functions in ArrayTraits.h.
 *
 * \todo Migrate to gsl_odeiv2 (although, it’s not clear what has changed there)
 */
template<typename A>
class MakerGSL<A>::_: public Evolved<A>
{
public:
  typedef Evolved<A> Base;

  typedef typename Base::Derivs Derivs;

  /** \see MakerGSL::MakerGSL() for an explanation of the role of `nextDtTryCorrectionFactor` */
  _(A&& a, Derivs derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs, SteppingFunction sf, double nextDtTryCorrectionFactor)
    : Base(std::forward<A>(a),derivs,dtInit,epsRel,epsAbs),
      pImpl_(details::createImpl(this,cpputils::size(this->getA()),auxFunction,epsRel,epsAbs,cpputils::data(scaleAbs),sf)),
      sf_(sf),
      nextDtTryCorrectionFactor_(nextDtTryCorrectionFactor)
  {
    if (!cpputils::isStorageContiguous(a)) throw (NonContiguousStorageException("In evolved"));
  }

  std::ostream& streamParameters_v(std::ostream& os) const {return os<<"EvolvedGSL implementation, stepping function: "<<sf_<<std::endl;}

private:

  static int auxFunction(double t, const double* y, double* dydt, void* aux)
  {
    using namespace cpputils;
    
    _* e=static_cast<_*>(aux);

    const A    yInterfaceA(create(y   ,e->getA()));
          A dydtInterfaceA(create(dydt,e->getA()));

    e->getDerivs()(t,yInterfaceA,dydtInterfaceA);

    return details::onSuccess;
  }

  void step_v(double deltaT, std::ostream&) final
  {
    double
      time=this->getTime(),
      dtTry=this->getDtTry(),
      nextDtTry=( fabs(deltaT)<fabs(dtTry/nextDtTryCorrectionFactor_) ? dtTry : 0. );

    apply(pImpl_,&time,time+deltaT,&dtTry,cpputils::data(this->getA()));

    Base::update(time, nextDtTry ? nextDtTry/nextDtTryCorrectionFactor_ : dtTry );

  }

  size_t nFailedStepsLast_v() const {return extractFailedSteps(pImpl_);}

  const details::ImplPtr pImpl_;

  const SteppingFunction sf_;
  const double nextDtTryCorrectionFactor_;

};




template<typename A>
auto MakerGSL<A>::make(
                       A&& a,
                       Derivs derivs,
                       double dtInit,
                       double epsRel,
                       double epsAbs,
                       const A& scaleAbs
                      ) const -> const Ptr
{
  return std::make_shared<_>(std::forward<A>(a),derivs,dtInit,epsRel,epsAbs,scaleAbs,sf_,nextDtTryCorrectionFactor_);
}



} // evolved

#endif // CPPQEDCORE_UTILS_EVOLVEDGSL_TCC_INCLUDED
