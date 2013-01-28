// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED
#define QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED

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

  typedef          structure::Averaged<RANK>      Averaged   ;
  typedef typename structure::Averaged<RANK>::Ptr AveragedPtr;

  DO_Display(
	     AveragedPtr,
	     const Pars&,
	     bool negativity,
	     size_t equalCount=10
	     ) throw(DimensionalityMismatchException);

  virtual ~DO_Display() {}

  std::ostream& displayMore   (double t, const DensityOperator&, std::ostream&, int precision) const throw(StoppingCriterionReachedException);
  size_t        displayMoreKey(std::ostream&) const;

private:
  const AveragedPtr av_ ;

  const bool negativity_;

  // Helpers for autoStopping
  const double autoStop_;

  mutable double lastCrit_;
  mutable size_t equalCount_;


};


} // details


} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_DO_DISPLAY_H_INCLUDED
