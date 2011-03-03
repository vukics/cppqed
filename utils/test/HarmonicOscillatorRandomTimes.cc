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
#include "Randomized.h"

#include <boost/bind.hpp>


using namespace std       ;
using namespace trajectory;
using namespace mathutils ;
using namespace parameters;


typedef TTD_CARRAY(1) CA1D;
typedef evolved::Evolved<CA1D>::SmartPtr ESP;

/*
  y(0) y
  y(1) dy/dt
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



void derivs(double tau, const CA1D& y, CA1D& dydt, DDOscPars& pars)
{
  pars.add();
  dydt(0)=y(1);
  dydt(1)=exp(DCOMP_I*pars.omega*tau)-2*pars.gamma*y(1)-y(0);
  // cout<<"# Derivs success"<<endl;
}

/*

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
  */  /*
  dcomp denom(sqr(p.omega)+2.*DCOMP_I*p.gamma*p.omega-1.)


  dcomp sqrtGamsqrMin1(sqrt(dcomp(p.gamma)-1.));
  dcomp denom(sqrtGamsqrMin1*());


  dcomp c1((y0*(1-p.gamma/sqrt));
      *//*

  return y0Block;

}
	*/


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsTrajectory pt(p);

  double 
    &omega=p.add("O","Driving frequency",1.),
    &gamma=p.add("G","Damping rate"     ,1.);

  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=10./min(1.,min(omega,gamma));

  double dtInit=1./(max(1.,max(omega,gamma))*TrajectoryBase::factor());  

  DDOscPars params(omega,gamma);

  CA1D y(2); y(0)=dcomp(1.,-1.); y(1)=dcomp(-1.,1.);

  Simulated<CA1D> simulated(y,bind(derivs,_1,_2,_3,params),dtInit,CA1D(),pt);
  simulated.displayParameters();

  randomized::Randomized::SmartPtr Ran(randomized::MakerGSL()(1001));

  for (size_t n=100; n>0; n--) {
    evolved::evolveTo(simulated,pt.T*(*Ran)()); simulated.display();
  }


  /*
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
  */
}
