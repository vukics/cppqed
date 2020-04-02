// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
// -*- C++ -*-
#ifndef CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED
#define CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED

#include "SimulatedFwd.h"

#include "ParsTrajectory.h"
#include "Trajectory.h"


namespace trajectory {

/// Class fully implementing the Adaptive interface by displaying (and serializing) the whole content of the evolved array
/**
 * Meant for all cases when simple ODE evolution is desired with intermittent displays
 * 
 * <b>Example usage:</b> simulation of a complex driven damped harmonic oscillator mode described by the ODE \f[\ddot{y}+2\gamma\,\dot{y}+y=e^{i\,\omega t},\f]
 * where \f$\gamma\f$ is the damping rate and \f$\omega\f$ the driving frequency, and the timescale has been chosen such that the eigenfrequency is 1.
 * 
 * \include HarmonicOscillatorComplex.cc
 * 
 * \todo Provide optional key printing
 */
template<typename A> 
class Simulated : public Adaptive<A>
{
public:
  typedef Adaptive<A> Base;

  typedef evolved::Evolved<A> Evolved;

  Simulated(A&, typename Evolved::Derivs, double dtInit,
            double, double,
            const A& scaleAbs=A(),
            const evolved::Maker<A>& =evolved::MakerGSL<A>());

  Simulated(A& array, typename Evolved::Derivs derivs, double dtInit,
            const ParsEvolved& pe,
            const A& scaleAbs=A(),
            const evolved::Maker<A>& maker=evolved::MakerGSL<A>()) : Simulated(array,derivs,dtInit,pe.epsRel,pe.epsAbs,scaleAbs,maker) {}

private:
  virtual void step_v(double deltaT) final {this->getEvolved()->step(deltaT);}

  virtual std::ostream& display_v(std::ostream&, int) const final;
  
  virtual std::ostream& displayKey_v(std::ostream& os, size_t&) const final {return os;}

  virtual std::ostream& displayParameters_v(std::ostream& os) const final {return Base::displayParameters_v(os<<"\n# Simulated.\n");}

  virtual const std::string trajectoryID_v() const final {return "Simulated";}

};


} // trajectory


#endif // CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED
