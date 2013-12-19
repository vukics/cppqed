// -*- C++ -*-
#ifndef   UTILS_SIMULATED_TCC_INCLUDED
#define   UTILS_SIMULATED_TCC_INCLUDED

#include "Simulated_.h"

// #include "ArrayTraits.h"
// The same note applies as with EvolvedGSL.tcc
#include "FormDouble.tcc"
#include "Trajectory.tcc"


namespace trajectory {


template<typename A>
Simulated<A>::Simulated(A& y, typename Evolved::Derivs derivs, double dtInit, 
                        double epsRel, double epsAbs,
                        const A& scaleAbs,
                        const evolved::Maker<A>& maker)
  : Base(y,derivs,dtInit,epsRel,epsAbs,scaleAbs,maker)
{}


template<typename A>
std::ostream& Simulated<A>::display_v(std::ostream& os, int precision) const
{
  using namespace cpputils;
  const A& a=getEvolved()->getA();
  for (size_t i=0; i<subscriptLimit(a); i++)
    os<<FormDouble(precision)(subscript(a,i))<<' ';
  return os;
}


} // trajectory

#endif // UTILS_SIMULATED_TCC_INCLUDED
