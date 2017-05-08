// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_EVOLVED_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_EVOLVED_TCC_INCLUDED

#include "Evolved.h"

#include "MathExtensions.h"

#include <boost/make_shared.hpp>

namespace evolved {

template<typename A>
EvolvedIO<A>::EvolvedIO(A& a, double dtInit, double epsRel, double epsAbs)
  : TimeStepBookkeeper(dtInit,epsRel,epsAbs), a_(a)
{}

template<typename A>
Evolved<A>::Evolved(A& a, Derivs derivs, double dtInit, double epsRel, double epsAbs) 
  : EvolvedIO<A>(a,dtInit,epsRel,epsAbs), derivs_(derivs)
{} 

template<typename A>
typename EvolvedIO<A>::Ptr makeIO(A& a, double time=0)
{
  typename EvolvedIO<A>::Ptr res = boost::make_shared<EvolvedIO<A>, A &>(a,0,0,0);
  res->setTime(time);
  return res;
}

template<typename A>
void Evolved<A>::step(double deltaT)
{
  if (mathutils::sign(deltaT)!=mathutils::sign(EvolvedIO<A>::getDtTry())) {
    // Stepping backward
    EvolvedIO<A>::setDtTry(-EvolvedIO<A>::getDtDid());
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

#endif // CPPQEDCORE_UTILS_EVOLVED_TCC_INCLUDED
