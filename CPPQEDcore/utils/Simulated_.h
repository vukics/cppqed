// -*- C++ -*-
#ifndef UTILS_SIMULATED__H_INCLUDED
#define UTILS_SIMULATED__H_INCLUDED

#include "SimulatedFwd.h"

#include "ParsTrajectory.h"
#include "Trajectory.h"


namespace trajectory {


/**
 * Example usage: \include HarmonicOscillatorComplex.cc
 * \todo Provide optional key printing
 */
template<typename A> 
class Simulated : public Adaptive<A>
{
public:
  typedef Adaptive<A> Base;

  typedef evolved::Evolved<A> Evolved;

  using Base::getEvolved;

  Simulated(A&, typename Evolved::Derivs, double dtInit,
            double, double,
            const A& scaleAbs=A(),
            const evolved::Maker<A>& =evolved::MakerGSL<A>());

  Simulated(A& array, typename Evolved::Derivs derivs, double dtInit,
            const ParsEvolved& pe,
            const A& scaleAbs=A(),
            const evolved::Maker<A>& maker=evolved::MakerGSL<A>()) : Simulated(array,derivs,dtInit,pe.epsRel,pe.epsAbs,scaleAbs,maker) {}

private:
  void step_v(double deltaT) {getEvolved()->step(deltaT);}

  std::ostream& display_v(std::ostream&, int) const;
  
  std::ostream& displayKey_v(std::ostream& os, size_t&) const {return os;}

  std::ostream& displayParameters_v(std::ostream& os) const {return Base::displayParameters_v(os<<"\n# Simulated.\n");}

  const std::string trajectoryID_v() const {return "Simulated";}

};


} // trajectory


#endif // UTILS_SIMULATED__H_INCLUDED
