// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolved.h"

#include <iostream>

namespace evolved {


TimeStepBookkeeper::TimeStepBookkeeper(double dtInit, double epsRel, double epsAbs)
  : t_(0), dtTry_(dtInit), dtDid_(0), epsRel_(epsRel), epsAbs_(epsAbs) 
{}

TimeStepBookkeeper::TimeStepBookkeeper(const TimeStepBookkeeper& t)
  : t_(t.t_), dtTry_(t.dtTry_), dtDid_(t.dtDid_), epsRel_(t.epsRel_), epsAbs_(t.epsAbs_) 
{}


void TimeStepBookkeeper::update(double t, double dtTry)
{
  dtDid_=t-t_ ;
  t_    =t    ;
  dtTry_=dtTry;
}


TimeStepBookkeeper& TimeStepBookkeeper::operator=(const TimeStepBookkeeper& other)
{
  if (&other!=this) {
    dtDid_=other.dtDid_;
    t_    =other.t_    ;
    dtTry_=other.dtTry_;
  }
  return *this;
}


using namespace std;

ostream& operator<<(ostream& os, SteppingFunction sf)
{
  switch (sf) {
  case SF_RKCK  : return os<<"rkck" ;
  case SF_RK8PD :        os<<"rk8pd";
  }
  return os;
}

istream& operator>>(istream& is, SteppingFunction& sf) 
{
  SteppingFunction sftemp=SF_RK8PD;
  string s;

  is>>s;
       if (s=="rkck"  ) sftemp=SF_RKCK ;
  else if (s!="rk8pd") 
    is.clear(ios_base::badbit);

  if (is) sf=sftemp;
  return is;
}


std::ostream& LoggingBase::logOnEnd(std::ostream& os) const
{
  return
  os<<"\nTotal number of ODE steps: "<<nSteps_
    <<"\nNumber of failed ODE steps: "<<nFailedSteps_
    <<"\nNumber of calls of function calculating RHS for ODE: "<<nDerivsCalls_<<std::endl;
}

} // evolved
