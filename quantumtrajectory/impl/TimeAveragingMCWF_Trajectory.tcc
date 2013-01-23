// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_TIMEAVERAGINGMCWF_TRAJECTORY_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_TIMEAVERAGINGMCWF_TRAJECTORY_TCC_INCLUDED

#include "TimeAveragingMCWF_Trajectory.h"

#include "impl/MCWF_Trajectory.tcc"


template<int RANK>
std::ostream& quantumtrajectory::TimeAveragingMCWF_Trajectory<RANK>::displayMore() const
{
  if (av_) {
    Averages averagesNow(av_->average(getTime(),getPsi()));
    if (getTime()>relaxationTime_) {
      averages_ = ( averages_*sum_ + averagesNow )/(sum_+1) ;
      ++sum_;
    }
    av_->process(averagesNow);
    av_->display(averagesNow,getOstream(),getPrecision());
  }
  return getOstream();
}


template<int RANK>
quantumtrajectory::TimeAveragingMCWF_Trajectory<RANK>::~TimeAveragingMCWF_Trajectory()
{
  if (av_) {
    std::ostream& os=getOstream()<<"# Time averages:\n# ";
    av_->process(averages_);
    av_->display(averages_,os,getPrecision())<<std::endl;
  }

}


#endif // QUANTUMTRAJECTORY_IMPL_TIMEAVERAGINGMCWF_TRAJECTORY_TCC_INCLUDED
