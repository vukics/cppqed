// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_TIMEAVERAGINGMCWF_TRAJECTORY_H_INCLUDED
#define   QUANTUMTRAJECTORY_TIMEAVERAGINGMCWF_TRAJECTORY_H_INCLUDED

#include "TimeAveragingMCWF_TrajectoryFwd.h"

#include "MCWF_Trajectory.h"


namespace quantumtrajectory {


template<int RANK>
class TimeAveragingMCWF_Trajectory : public MCWF_Trajectory<RANK>
{
public:
  typedef MCWF_Trajectory<RANK> Base;
  
  typedef typename Base::StateVector    StateVector   ;
  typedef typename Base::StateVectorLow StateVectorLow;
  
  typedef typename Base::Averaged Averaged;
  
  typedef typename Averaged::Averages Averages;
  
  using Base::getQS; using Base::getDtDid; using Base::getTime; using Base::getPsi;
  
  template<typename SYS>
  TimeAveragingMCWF_Trajectory(
                               StateVector& psi,
                               const SYS& sys,
                               const ParsMCWF& p,
                               double relaxationTime,
                               const StateVectorLow& scaleAbs=StateVectorLow()
                               )
    : MCWF_Trajectory<RANK>(psi,sys,p,scaleAbs), relaxationTime_(relaxationTime), averages_(getQS().template nAvr<structure::LA_Av>()), sum_(0), av_(getQS().getAv())
    {
      averages_=0.;
    }

  ~TimeAveragingMCWF_Trajectory();

  const Averages getAverages() const {return averages_;}
  
private:
  std::ostream& display_v(std::ostream&, int) const;
  
  const double relaxationTime_;
  
  mutable Averages averages_;
  mutable double sum_;
  
  const typename Averaged::Ptr av_;

};


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_TIMEAVERAGINGMCWF_TRAJECTORY_H_INCLUDED
