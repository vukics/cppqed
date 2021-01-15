// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED
#define CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED

#include "ParsTrajectory.h"
#include "Trajectory.h"
// #include "ArrayTraits.h"
// The same note applies as with EvolvedGSL.tcc
#include "FormDouble.h"
#include "Trajectory.tcc"


namespace trajectory {

/// Class fully implementing the Adaptive interface by streaming (and serializing) the whole content of the evolved array
/**
 * Meant for all cases when simple ODE evolution is needed with intermittent streamings
 * 
 * <b>Example usage:</b> simulation of a complex driven damped harmonic oscillator mode described by the ODE \f[\ddot{y}+2\gamma\,\dot{y}+y=e^{i\,\omega t},\f]
 * where \f$\gamma\f$ is the damping rate and \f$\omega\f$ the driving frequency, and the timescale has been chosen such that the eigenfrequency is 1.
 * 
 * \include HarmonicOscillatorComplex.cc
 * 
 * \todo Provide optional key printing
 */
template<typename A> 
class Simulated : public Adaptive<A,Trajectory<A> >
{
public:
  Simulated(Simulated&&) = default; Simulated& operator=(Simulated&&) = default;
  
  typedef Adaptive<A,Trajectory<A> > Base;

  typedef evolved::Evolved<A> Evolved;

  template <typename ARRAY>
  Simulated(ARRAY&& y, typename Evolved::Derivs derivs, double dtInit,
            int logLevel,
            double epsRel, double epsAbs,
            const A& scaleAbs=A(),
            const evolved::Maker<A>& maker=evolved::MakerGSL<A>())
    : Base(std::forward<ARRAY>(y),derivs,dtInit,logLevel,epsRel,epsAbs,scaleAbs,maker) {}

  template <typename ARRAY>
  Simulated(ARRAY&& array, typename Evolved::Derivs derivs, double dtInit,
            const ParsEvolved& pe,
            const A& scaleAbs=A(),
            const evolved::Maker<A>& maker=evolved::MakerGSL<A>()) : Simulated(std::forward<ARRAY>(array),derivs,dtInit,pe.logLevel,pe.epsRel,pe.epsAbs,scaleAbs,maker) {}

protected:
  typename Base::StreamReturnType stream_v(std::ostream& os, int precision) const override
  {
    using namespace cpputils;
    const A& a=this->getEvolved()->getA();
    for (size_t i=0; i<subscriptLimit(a); i++) os<<FormDouble(precision)(subscript(a,i))<<' ';
    return {os,a.copy()};
  }
  
private:
  void step_v(double deltaT, std::ostream& logStream) final {this->getEvolved()->step(deltaT,logStream);}

  std::ostream& streamKey_v(std::ostream& os, size_t&) const final {return os;}

  std::ostream& streamParameters_v(std::ostream& os) const final {return Base::streamParameters_v(os<<"\nSimulated.\n");}

  const std::string trajectoryID_v() const final {return "Simulated";}

};


} // trajectory


#endif // CPPQEDCORE_UTILS_SIMULATED__H_INCLUDED
