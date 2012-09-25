#include "Simulated.h"

#include "MathExtensions.h"
#include "Randomized.h"

#include <boost/bind.hpp>


using namespace std       ;
using namespace trajectory;
using namespace mathutils ;
using namespace parameters;


typedef TTD_CARRAY(1) CA1D;
typedef evolved::Evolved<CA1D>::Ptr ESP;

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
}


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsTrajectory pt(p);

  double 
    &omega=p.add("O","Driving frequency",1.),
    &gamma=p.add("G","Damping rate"     ,1.);

  pt.epsRel=1e-12;

  // Parameter finalization
  update(p,argc,argv,"--");

  if (pt.T<0) pt.T=20./min(1.,min(omega,gamma));

  double dtInit=1./(max(1.,max(omega,gamma))*TrajectoryBase::factor());  

  DDOscPars params(omega,gamma);

  CA1D y(2); y(0)=dcomp(3.45263,-2.0746); y(1)=dcomp(-4.83065,1.16527);

  Simulated<CA1D> simulated(y,bind(derivs,_1,_2,_3,params),dtInit,CA1D(),pt);
  simulated.displayParameters();

  randomized::Randomized::Ptr Ran(randomized::MakerGSL()(1001));

  for (size_t n=1000; n>0; n--) {
    evolved::evolveTo(simulated,pt.T*(*Ran)()); simulated.display();
  }


}
