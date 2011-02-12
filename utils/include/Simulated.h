// -*- C++ -*-
#ifndef _SIMULATED_H
#define _SIMULATED_H

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

  void   displayMore(int) const;
  size_t displayMoreKey() const;

  void   displayParameters() const {getOstream()<<std::endl<<"# Simulated."<<std::endl; Base::displayParameters();}

};


} // trajectory

#include "impl/Simulated.tcc"

#endif // _SIMULATED_H
