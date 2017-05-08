// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_SIMULATED_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_SIMULATED_TCC_INCLUDED

#include "Simulated_.h"

// #include "ArrayTraits.h"
// The same note applies as with EvolvedGSL.tcc
#include "FormDouble.tcc"
#include "Trajectory.tcc"


namespace trajectory {


template<typename A>
Simulated<A>::Simulated(A& y, typename Evolved::Derivs derivs, double dtInit,
                        int logLevel,
                        double epsRel, double epsAbs,
                        const A& scaleAbs,
                        const evolved::Maker<A>& maker)
  : Base(y,derivs,dtInit,logLevel,epsRel,epsAbs,scaleAbs,maker)
{}


template<typename A>
std::ostream& Simulated<A>::display_v(std::ostream& os, int precision) const
{
  displayPreHook();
  using namespace cpputils;
  const A& a=this->getEvolved()->getA();
  for (size_t i=0; i<subscriptLimit(a); i++)
    os<<FormDouble(precision)(subscript(a,i))<<' ';
  displayPostHook();
  return os;
}


} // trajectory

#endif // CPPQEDCORE_UTILS_SIMULATED_TCC_INCLUDED
