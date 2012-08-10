/*
  
  A rather artificial example based on the observation that

  exp(exp(Zt))

  is a particular solution of the equation

  D^2(y) - D(y)^2 - Z D(y) + Z^2 exp(2Zt) ( y^2 - y ) = 0

  for initial condition

  y(0) = e
  p(0) = Z e

 */

#include "Simulated.h"

#include "MathExtensions.h"

#include<boost/bind.hpp>


using namespace std;
// using namespace blitzplusplus;
using namespace mathutils;
using namespace parameters;
using namespace trajectory;
/*
  y[0] Re y
  y[1] Im y
  y[2] Re dy/dt
  y[3] Im dy/dt
*/

// #include<cmath>
// #include<iostream>

typedef TTD_CARRAY(1) CA1D;
typedef TTD_DARRAY(1) DA1D;

void derivs(double tau, const CA1D& yA, CA1D& dydtA, const dcomp& Z)
{
  // static int count=0;

  dcomp
    y(yA(0)),
    p(yA(1));

  dydtA(0)=p;

  dydtA(1)=sqr(p)+Z*p+sqr(Z)*exp(2.*Z*tau)*(y-sqr(y));

  // cout<<count++<<endl;

}


int main()
{

  dcomp Z(-1.,10);

  stringstream slate;

  {
    double 
      T     =4.,//10./min(fabs(real(Z)),fabs(imag(Z))),
      dtInit=1./(abs(Z)*TrajectoryBase::factor());
 
    CA1D y(2);
  
    y(0)=  exp(1);
    y(1)=Z*exp(1);

    DA1D yr(real(y));

    Simulated<CA1D> S(y,bind(derivs,_1,_2,_3,Z),dtInit,1e-6,1e-30,CA1D(),slate,6);
  
    run(S,T,1,false);
  }

  double diff=0, sum=0;

  while (slate) {
    string str; getline(slate,str);
    if (str!="" && str[0]!='#') {
      stringstream line(str);
      double time;
      line>>time;
      dcomp 
	calculated(exp(exp(Z*time))),
	read;
      line>>time>>read; // eat the timestep first
      diff+=abs(read-calculated); sum+=(abs(read)+abs(calculated))/2.;
    }
  }

  return diff/sum<1e-5 ? 0 : 1;
  
}
