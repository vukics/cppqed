// -*- C++ -*-
#ifndef   _EVOLVED_IMPL_H
#define   _EVOLVED_IMPL_H


namespace evolved {


template<typename A>
Evolved<A>::Evolved(A& a, Derivs derivs, double dtInit, double epsRel, double epsAbs) 
  : EvolvedCommon(dtInit,epsRel,epsAbs), a_(a), derivs_(derivs)
{} 


template<typename E>
void evolve(E& e, double deltaT)
{
  double endTime=e.getTime()+deltaT;
  while (double dt=endTime-e.getTime()) e.step(dt);
}
// evolves for exactly deltaT


}

#endif // _EVOLVED_IMPL_H
