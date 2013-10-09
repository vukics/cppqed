#include "Evolved.h"

namespace evolved {


TimeStepBookkeeper::TimeStepBookkeeper(double dtInit, double epsRel, double epsAbs)
  : t_(0), dtTry_(dtInit), dtDid_(0), epsRel_(epsRel), epsAbs_(epsAbs) 
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



} // evolved
