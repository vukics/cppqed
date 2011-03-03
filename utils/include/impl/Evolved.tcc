// -*- C++ -*-
#ifndef   _EVOLVED_IMPL_H
#define   _EVOLVED_IMPL_H

#include "MathExtensions.h"

namespace evolved {


template<typename A>
Evolved<A>::Evolved(A& a, Derivs derivs, double dtInit, double epsRel, double epsAbs) 
  : EvolvedCommon(dtInit,epsRel,epsAbs), a_(a), derivs_(derivs)
{} 


template<typename A>
void Evolved<A>::step(double deltaT)
{
  if (mathutils::sign(deltaT)!=mathutils::sign(getDtTry())) {
    // Stepping backward
    setDtTry(-getDtDid());
    doStep(deltaT);
  }
  else doStep(deltaT);
}


template<typename E>
void evolve(E& e, double deltaT)
{
  double endTime=e.getTime()+deltaT;
  while (double dt=endTime-e.getTime()) e.step(dt);
}
// evolves for exactly deltaT


}

#endif // _EVOLVED_IMPL_H
