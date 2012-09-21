#include "Simulated.h"

#include <boost/bind.hpp>


using namespace std       ;
using namespace trajectory;
using namespace parameters;


typedef blitz::Array<dcomp,1> CA1R;

/*
  y(0) y
  y(1) dy/dt
*/

void derivs(double tau, const CA1R& y, CA1R& dydt, double omega, double gamma)
{
  dydt(0)=y(1);
  dydt(1)=exp(DCOMP_I*omega*tau)-2*gamma*y(1)-y(0);
}


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsTrajectory pt(p);

  double& omega=p.add("O","Driving frequency",1.);
  double& gamma=p.add("G","Damping rate"     ,1.);

  dcomp&    yinit=p.add(   "yinit"," y   initial condition",dcomp( 3.45263,-2.0746 ));
  dcomp& dydtinit=p.add("dydtinit","dydt initial condition",dcomp(-4.83065, 1.16527));

  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=20./min(1.,min(omega,gamma));
  // Note: 1.0 is also and existing frequency in the system, which defines the unit of time 

  CA1R y(2); y=yinit,dydtinit;

  Simulated<CA1R> S(y,
		    bind(derivs,_1,_2,_3,omega,gamma),
		    .1/max(1.,max(omega,gamma)),
		    CA1R(),
		    pt);

  evolve(S,pt);

}
