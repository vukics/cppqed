// -*- C++ -*-
#ifndef _DO_DISPLAY_H
#define _DO_DISPLAY_H

#include "DensityOperatorFwd.h"
#include "DO_DisplayFwd.h"

#include "DimensionsBookkeeper.h"

#include "AveragedFwd.h"
#include "QuantumSystemFwd.h"

#include "Trajectory.h"

// #include "NegPT.h"


namespace quantumtrajectory {


namespace details {


using namespace trajectory;


template<int RANK, typename V>
class DO_Display 
{
public:
  typedef quantumdata::DensityOperator<RANK> DensityOperator;

  typedef structure::QuantumSystem<RANK> QuantumSystem;
  typedef structure::Averaged     <RANK> Averaged     ;  

  DO_Display(
	     const QuantumSystem&, 
	     const ParsTrajectory&,
	     bool negativity,
	     size_t equalCount=10
	     ) throw(DimensionalityMismatchException);

  virtual ~DO_Display() {}

  const Averaged*const getAveraged() const {return av_;}

  void   displayMore(double t, const DensityOperator&, std::ostream&, int precision) const throw(StoppingCriterionReachedException);
  size_t displayMoreKey(std::ostream&) const;

private:
  const Averaged*const av_ ;

  const bool negativity_;

  // Helpers for autoStopping
  const double autoStop_;

  mutable double lastCrit_;
  mutable size_t equalCount_;


};


} // details


} // quantumtrajectory

#include "impl/DO_Display.tcc"

#endif // _DO_DISPLAY_H
