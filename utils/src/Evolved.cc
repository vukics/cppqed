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


} // evolved
