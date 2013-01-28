// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_H_INCLUDED
#define   QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_H_INCLUDED

#include "ProjectingMCWF_TrajectoryFwd.h"

#include "MCWF_Trajectory.h"

#include <boost/ptr_container/ptr_vector.hpp>


namespace quantumtrajectory {


template<int RANK>
class ProjectingMCWF_Trajectory : public MCWF_Trajectory<RANK>
{
public:
  typedef MCWF_Trajectory<RANK> Base;

  typedef typename Base::StateVector    StateVector   ;
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef boost::ptr_vector<StateVector> Basis;

  using Base::getOstream; using Base::getPsi; using Base::getPrecision;

  template<typename SYS>
  ProjectingMCWF_Trajectory(
			    StateVector& psi,
			    const Basis& basis,
			    const SYS& sys,
			    const ParsMCWF_Trajectory& p,
			    const StateVectorLow& scaleAbs=StateVectorLow()
			    )
    : trajectory::Trajectory(p), Base(psi,sys,p,scaleAbs), basis_(basis), metricTensor_uu_(help())
  {}

private:
  std::ostream& displayMore() const;
  size_t displayMoreKey() const;

  const linalg::CMatrix help() const;

  const Basis basis_;
  const linalg::CMatrix metricTensor_uu_;

};


} // quantumtrajectory



#endif // QUANTUMTRAJECTORY_PROJECTINGMCWF_TRAJECTORY_H_INCLUDED
