// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_TIMEAVERAGINGMCWF_TRAJECTORY_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_TIMEAVERAGINGMCWF_TRAJECTORY_TCC_INCLUDED

#include "TimeAveragingMCWF_Trajectory.h"

#include "impl/MCWF_Trajectory.tcc"


template<int RANK>
std::ostream& quantumtrajectory::TimeAveragingMCWF_Trajectory<RANK>::displayMore() const
{
  if (getTime()>relaxationTime_) {
    const double dtDid=getDtDid();
    averages_ = ( averages_*sumDt_ + getQS().template average<structure::LA_Av>(getTime(),getPsi())*dtDid )/(sumDt_+dtDid) ;
    sumDt_+=dtDid;
  }
  return Base::displayMore();
}


template<int RANK>
quantumtrajectory::TimeAveragingMCWF_Trajectory<RANK>::~TimeAveragingMCWF_Trajectory()
{
  const typename Averaged::Ptr av(getQS().getAv());
  if (av) {
    std::ostream& os=getOstream()<<"# Time averages:\n# ";
    av->process(averages_);
    av->display(averages_,os,getPrecision())<<std::endl;
  }

}


#endif // QUANTUMTRAJECTORY_IMPL_TIMEAVERAGINGMCWF_TRAJECTORY_TCC_INCLUDED
