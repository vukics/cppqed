// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_EVOLVEDNR_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_EVOLVEDNR_TCC_INCLUDED

#include "EvolvedNR.h"

namespace evolved {


template<typename A>
const typename Maker<A>::SmartPtr 
MakerNR<A>::operator()(A& a,
		       Derivs derivs,
		       double dtInit,
		       double epsRel,
		       double epsAbs, 
		       const A& scaleAbs
		       ) const
{
  return SmartPtr(new details::EvolvedNR<A>(a,derivs,dtInit,epsRel,epsAbs,scaleAbs));
}


namespace details {


template<typename A>
EvolvedNR<A>::EvolvedNR(
			A& a,
			Derivs derivs,
			double dtInit,
			double epsRel,
			double epsAbs,
			const A& scaleAbs
			)
  : Base(a,derivs,dtInit,epsRel,epsAbs/epsRel),
    yscal_(Traits::create(a)), dydt_(Traits::create(a)), scaleAbs_(scaleAbs)
{
}


template<typename A>
void 
rkqs(A&, const A&, double&, double, double, const A&, double&, double&, typename Evolved<A>::Derivs);


template<typename A>
void
EvolvedNR<A>::doStep(double deltaT)
{
  double t=getTime(), dtTry=getDtTry(), dtDummy;

  Base::getDerivs()(t,getA(),dydt_);

  if (scaleAbs_.size())
    yscal_=abs(getA())+abs(dtTry*dydt_)+Base::getEpsAbs()*scaleAbs_;
  else
    yscal_=abs(getA())+abs(dtTry*dydt_);

  rkqs(getA(),
       dydt_,
       t,
       std::min(dtTry,deltaT),
       Base::getEpsRel(),
       yscal_,
       dtDummy,
       dtTry,
       Base::getDerivs());
  Base::update(t,dtTry);

}


template<typename A>
void
rkqs(A& y, const A& dydx, double& x, double htry, double eps, const A& yscal, double& hdid, double& hnext,
     typename Evolved<A>::Derivs derivs)
{
  typedef cpputils::ArrayMemoryTraits<A> Traits;

  const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;

  double h=htry, errmax;

  A yerr(Traits::create(y)), ytemp(Traits::create(y));

  for (;;) {

    rkck(y,dydx,x,h,ytemp,yerr,derivs);

    errmax=max(abs(yerr/yscal))/eps;

    if (errmax <= 1.0) break;

    double htemp=SAFETY*h*pow(errmax,PSHRNK);
    h=(h >= 0.0 ? std::max(htemp,0.1*h) : std::min(htemp,0.1*h));

  }

  if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
  else hnext=5.0*h;
  x += (hdid=h);

  y=ytemp;

}

namespace cashkarpconsts {

const double 
a2=.2,
  a3=.3, 
  a4=.6, 
  a5=1., 
  a6=.875;

const double 
b21=.2;

const double
b31=3./40.,
  b32=9./40.;

const double
b41=.3,
  b42=-.9,
  b43=1.2;

const double
b51=-11./54.,
  b52=2.5,
  b53=-70./27.,
  b54=35./27.;

const double
b61=1631./55296.,
  b62=175./512.,
  b63=575./13824.,
  b64=44275./110592., 
  b65=253./4096.;
  
const double
c1=37./378.,
  c3=250./621.,
  c4=125./594.,
  c6=512./1771.;

const double
dc1=c1-2825./27648.,
  dc3=c3-18575./48384.,
  dc4=c4-13525./55296., 
  dc5=-277./14336.,
  dc6=c6-.25;

} // cashkarpconsts


template<typename A>
void 
rkck(const A& y, const A& dydx, double x, double h, A& yout, A& yerr, typename Evolved<A>::Derivs derivs)
{
  typedef cpputils::ArrayMemoryTraits<A> Traits;
  using namespace cashkarpconsts;

  A
    ak2(Traits::create(y)),
    ak3(Traits::create(y)),
    ak4(Traits::create(y)),
    ak5(Traits::create(y)),
    ak6(Traits::create(y)),
    ytemp(Traits::create(y));

  ytemp=y+b21*h*dydx;

  derivs(x+a2*h,ytemp,ak2);

  ytemp=y+h*(b31*dydx+b32*ak2);

  derivs(x+a3*h,ytemp,ak3);

  ytemp=y+h*(b41*dydx+b42*ak2+b43*ak3);

  derivs(x+a4*h,ytemp,ak4);

  ytemp=y+h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4);

  derivs(x+a5*h,ytemp,ak5);

  ytemp=y+h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5);

  derivs(x+a6*h,ytemp,ak6);

  yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6);

  yerr=h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);

}


} // details



} // evolved


#endif // UTILS_INCLUDE_IMPL_EVOLVEDNR_TCC_INCLUDED
