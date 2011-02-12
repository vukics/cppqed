/*
  
  Solution of the ODE

  diff(y(t),t,t)+2*gamma*diff(y(t),t)+y(t)-exp(-I*omega*t)=0

  gamma==1 (critically damped)
  
  y(t)=exp(-t)*_C2+exp(-t)*t*_C1-exp(-I*omega*t)/(-1+(2*I)*omega+omega^2)

  with
  _C2=(-y0+(2*I)*y0*omega+y0*omega^2+1)/(-1+(2*I)*omega+omega^2)
  _C1=((2*I)*Dy0*omega+(2*I)*y0*omega-I*omega+Dy0*omega^2+y0*omega^2-Dy0+1-y0)/(-1+(2*I)*omega+omega^2)

  gamma!=1

  y(t)=exp((-gamma+sqrt(gamma^2-1))*t)*_C2+exp((-gamma-sqrt(gamma^2-1))*t)*_C1-exp(-I*omega*t)/(-1+(2*I)*gamma*omega+omega^2)  

 */

#include "SimulatedHigh.h"

#include "MathExtensions.h"

#include<boost/bind.hpp>


using namespace std       ;
using namespace trajectory;
using namespace mathutils ;
using namespace parameters;


typedef TTD_DARRAY(1) DA1D;

/*
  y(0) Re y
  y(1) Im y
  y(2) Re dy/dt
  y(3) Im dy/dt
*/

class DDOscPars {

  size_t count_;

public:
  double omega;
  double gamma;

  DDOscPars(double o, double g) : count_(0), omega(o), gamma(g) {}

  void reset() {count_=0;}
  void add() {count_++;}
  size_t count() const {return count_;}

};



void derivs(double tau, const DA1D& y, DA1D& dydt, DDOscPars& pars)
{
  pars.add();
  dydt(0)=y(2); dydt(1)=y(3); 
  dydt(2)=cos(pars.omega*tau)-2*pars.gamma*y(2)-y(0); dydt(3)=sin(pars.omega*tau)-2*pars.gamma*y(3)-y(1);
  cout<<"# Derivs success"<<endl;
}

DA1D& initial(DA1D& y)
{
  y(0)=1.; y(1)=-1.; y(2)=-1.; y(3)=1.;
  return y;
}

DA1D& exactSolution(const DDOscPars& p, double t, DA1D& y0Block)
{
  dcomp y0
    (y0Block(0),y0Block(1));
  dcomp Dy0
    (y0Block(2),y0Block(3));

  if (p.gamma==1.) {

    dcomp denom((-1.+(2.*DCOMP_I)*p.omega+sqr(p.omega)));

    dcomp 
      _C2(
	  (-y0+(2.*DCOMP_I)*y0*p.omega+y0*sqr(p.omega)+1.)
	  /
	  denom),
      _C1(((2.*DCOMP_I)*Dy0*p.omega+(2.*DCOMP_I)*y0*p.omega-DCOMP_I*p.omega+Dy0*sqr(p.omega)+y0*sqr(p.omega)-Dy0+1.-y0)
	  /
	  denom);

    // cout<<_C2<<' '<<_C1<<endl;
    

    dcomp 
      y(exp(-t)*_C2+exp(-t)*t*_C1-exp(-DCOMP_I*p.omega*t)
	/
	denom);
    
    y0Block(0)=real(y);
    y0Block(1)=imag(y);
    
  }
  /*
  dcomp denom(sqr(p.omega)+2.*DCOMP_I*p.gamma*p.omega-1.)


  dcomp sqrtGamsqrMin1(sqrt(dcomp(p.gamma)-1.));
  dcomp denom(sqrtGamsqrMin1*());


  dcomp c1((y0*(1-p.gamma/sqrt));
  */

  return y0Block;

}

/*

Evolved& EM(int mode, DA1D& y, double dtInit, void* params, double eps)
{
  switch (mode) {
  default: 
    return GSL_EvolvedMaker().make(y,4,derivs,dtInit,params,eps);
  case 1:
    return *(new GSL_Evolved(y,4,derivs,dtInit,params,eps,gsl_odeiv_step_rk8pd));
  }
}

*/


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsTrajectory pt(p);

  double& omega=p.add("O","Driving frequency",1.);
  double& gamma=p.add("G","Damping rate"     ,1.);

  bool &exact=p.add("exact","Exact solution switch",false);
  
  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=100./min(1.,min(omega,gamma));
  double dtInit=.1/(max(1.,max(omega,gamma))*TrajectoryBase::factor());  

  DDOscPars params(omega,gamma);

  DA1D y(4); 

  if (exact) {
    for (double t=0; t<pt.T; t+=dtInit) {
      y=exactSolution(params,t,initial(y));
      cout<<t<<' '<<0<<' '<<y(0)<<' '<<y(1)<<endl;
    }
  }
  else {
    initial(y);

    Simulated<DA1D> S(y,bind(derivs,_1,_2,_3,params),dtInit,DA1D(),pt);

    if (pt.T<0) pt.T=100./min(1.,min(omega,gamma));

    run(S,pt.T,pt.dc,true);

  }

}
