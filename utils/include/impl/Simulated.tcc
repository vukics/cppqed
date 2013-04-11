// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_SIMULATED_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_SIMULATED_TCC_INCLUDED

#include "Simulated_.h"

// #include "ArrayTraits.h"
// The same note applies as with EvolvedGSL.tcc
#include "impl/FormDouble.tcc"
#include "impl/Trajectory.tcc"


namespace trajectory {


template<typename A>
Simulated<A>::Simulated(A& y, typename Evolved::Derivs derivs, double dtInit, 
                        double epsRel, double epsAbs,
                        const A& scaleAbs, 
                        std::ostream& os, int precision, 
                        const evolved::Maker<A>& maker)
  : Base(y,derivs,dtInit,epsRel,epsAbs,scaleAbs,maker)
{}


template<typename A> 
Simulated<A>::Simulated(A& y, typename Evolved::Derivs derivs, double dtInit,
                        const A& scaleAbs,
                        const ParsEvolved& p,
                        const evolved::Maker<A>& maker)
  : Base(y,derivs,dtInit,p,scaleAbs,maker)
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

#endif // UTILS_INCLUDE_IMPL_SIMULATED_TCC_INCLUDED
