// -*- C++ -*-
#ifndef   UTILS_EVOLVED_TCC_INCLUDED
#define   UTILS_EVOLVED_TCC_INCLUDED

#include "Evolved.h"

#include "MathExtensions.h"


namespace evolved {


template<typename A>
Evolved<A>::Evolved(A& a, Derivs derivs, double dtInit, double epsRel, double epsAbs) 
  : TimeStepBookkeeper(dtInit,epsRel,epsAbs), a_(a), derivs_(derivs)
{} 


template<typename A>
void Evolved<A>::step(double deltaT)
{
  if (mathutils::sign(deltaT)!=mathutils::sign(getDtTry())) {
    // Stepping backward
    setDtTry(-getDtDid());
    step_v(deltaT);
  }
  else step_v(deltaT);
}



template<typename E>
void evolve(E& e, double deltaT)
{
  double endTime=e.getTime()+deltaT;
  while (double dt=endTime-e.getTime()) e.step(dt);
}
// evolves for exactly deltaT


} // evolved

#endif // UTILS_EVOLVED_TCC_INCLUDED
