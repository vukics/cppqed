// -*- C++ -*-
#ifndef UTILS_INCLUDE_SIMULATED__H_INCLUDED
#define UTILS_INCLUDE_SIMULATED__H_INCLUDED

#include "SimulatedFwd.h"

#include "ArrayTraitsFwd.h"
#include "EvolvedFwd.h"

#include "Trajectory.h"


namespace trajectory {


template<typename A> 
class Simulated : public Trajectory<A>
{
public:
  typedef Trajectory<A> Base;
  typedef cpputils::ArrayTraversalTraits<A> Traits;

  typedef evolved::Evolved<A> Evolved;

  using Base::getEvolved; using Base::getOstream;

  Simulated(A&, typename Evolved::Derivs, double dtInit, 
	    double, double, 
	    const A& scaleAbs,
	    std::ostream&, int,
	    const evolved::Maker<A>& =evolved::MakerGSL<A>());

  Simulated(A&, typename Evolved::Derivs, double dtInit,
	    const A& scaleAbs,
	    const ParsTrajectory&,
	    const evolved::Maker<A>& =evolved::MakerGSL<A>());

  void step(double deltaT) const {getEvolved()->step(deltaT);}

  void   displayMore   () const;
  size_t displayMoreKey() const;

  void   displayParameters() const {getOstream()<<std::endl<<"# Simulated."<<std::endl; Base::displayParameters();}

};


} // trajectory


#endif // UTILS_INCLUDE_SIMULATED__H_INCLUDED
