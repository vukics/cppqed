// -*- C++ -*-
#ifndef UTILS_INCLUDE_SIMULATED__H_INCLUDED
#define UTILS_INCLUDE_SIMULATED__H_INCLUDED

#include "SimulatedFwd.h"

#include "EvolvedFwd.h"

#include "Trajectory.h"


namespace trajectory {


template<typename A> 
class Simulated : public Adaptive<A>
{
public:
  typedef Adaptive<A> Base;

  typedef evolved::Evolved<A> Evolved;

  using Base::getEvolved;

  Simulated(A&, typename Evolved::Derivs, double dtInit,
            double, double,
            const A& scaleAbs,
            const evolved::Maker<A>& =evolved::MakerGSL<A>());

  Simulated(A&, typename Evolved::Derivs, double dtInit,
            const A& scaleAbs,
            const ParsEvolved&,
            const evolved::Maker<A>& =evolved::MakerGSL<A>());

private:
  void step_v(double deltaT) {getEvolved()->step(deltaT);}

  std::ostream& display_v(std::ostream&, int) const;
  
  std::ostream& displayKey_v(std::ostream& os, size_t&) const {return os;}

  std::ostream& displayParameters_v(std::ostream& os) const {return Base::displayParameters_v(os<<"\n# Simulated.\n");}

  std::string trajectoryID_v() const {return trajectoryID_;}
  static const char trajectoryID_[];
};


} // trajectory


#endif // UTILS_INCLUDE_SIMULATED__H_INCLUDED
