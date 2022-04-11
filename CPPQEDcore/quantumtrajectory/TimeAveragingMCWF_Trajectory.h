// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_TIMEAVERAGINGMCWF_TRAJECTORY_H_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_TIMEAVERAGINGMCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_Trajectory.h"


namespace quantumtrajectory {


/// Derived from MCWF_Trajectory, this class in addition performs time averaging of the columns during trajectory stream
/**
 * The time averages and the amount of averaging are part of the trajectory state, so that time averaging can be correctly resumed over trajectory-state i/o.
 * Averaging starts after a specified relaxation time.
 * 
 * \note Time averaging does not use stepsize-weighting, as experience has shown that this leads to worse convergence (similarly to trajectory::Ensemble).
 * 
 * \todo Stepsize-weighting could eventually be enabled as an option by a switch
 * 
 * \par An important note concerning sampling
 * Let us identify the two opposed regimes of our adaptive MCWF evolution:
 * 1. **weak noise**: when the contribution of jumps is small compared to the coherent development, the timestep-control is ODE-dominated, we will hardly ever experience dpLimit-overshoots
 * 2. **strong noise**: in the opposite regime, the timestep-control is dpLimit-dominated, the timestep will be approximately dpLimit/(total jump rate)
 * 
 * Let us consider a harmonic-oscillator mode, where the jump rate is a monotonic function of the photon number. In the 2. case, there will hence be perfect correlation between the 
 * photon number and the stepsize, which will lead to erroneous sampling in dc-mode, when the samples with the smaller stepsize are represented with the larger weight in the stream.
 * 
 * Hence, in this case one must be careful to use deltaT mode.
 * 
 * \see The different versions of trajectory::run
 * 
 * \todo If need be, the class could be modified to accept also a ProjectingMCWF_Trajectory as a base
 */
template<int RANK>
class TimeAveragingMCWF_Trajectory : public MCWF_Trajectory<RANK>
{
private:
  typedef MCWF_Trajectory<RANK> Base;
  
public:
  typedef typename Base::StateVector    StateVector   ;
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::SV_Ptr SV_Ptr;
  
  typedef typename Base::Averaged Averaged;
  
  typedef typename Averaged::Averages Averages;

  /// The signature is identical to MCWF_Trajectory::MCWF_Trajectory, but the relaxation time must be supplied as well.
  TimeAveragingMCWF_Trajectory(
                               SV_Ptr psi,
                               typename structure::QuantumSystem<RANK>::Ptr sys,
                               const mcwf::Pars& p,
                               double relaxationTime, ///< relaxation time after which the time averaging starts
                               const StateVectorLow& scaleAbs=StateVectorLow()
                               )
    : Base(psi,sys,p,scaleAbs), relaxationTime_(relaxationTime), averages_(this->template nAvr<structure::LA_Av>()), sum_(0), av_(this->getAv())
    {
      averages_=0.;
    }

  const Averages getAverages() const {return averages_;}
  
private:
  typename Base::StreamReturnType stream_v(std::ostream& os, int precision) const override
  {
    Averages averagesNow;
    if (av_) {
      averagesNow.reference(av_->average(this->getTime(),*this->getPsi()));
      if (this->getTime()>relaxationTime_) {
        averages_ = ( averages_*sum_ + averagesNow )/(sum_+1) ;
        ++sum_;
      }
      av_->process(averagesNow);
      av_->stream(averagesNow,os,precision);
    }
    return {os,averagesNow};
  }
  
  std::ostream& logMoreOnEnd_v(std::ostream& os) const override
  {
    if (av_) {
      Averages endValues(averages_.copy());
      // process would “corrupt” the averages_, which would be incorrectly saved afterwards by writeState_v, therefore we process a copy instead
      av_->process(endValues);
      av_->stream(endValues,os<<"Time averages after relaxation time "<<relaxationTime_<<std::endl,FormDouble::overallPrecision)<<std::endl;
    }
    return Base::logMoreOnEnd_v(os);
  }
  
  cppqedutils::iarchive&  readStateMore_v(cppqedutils::iarchive& iar)       override {return Base:: readStateMore_v(iar) & averages_ & sum_;}
  cppqedutils::oarchive& writeStateMore_v(cppqedutils::oarchive& oar) const override {return Base::writeStateMore_v(oar) & averages_ & sum_;}
  
  const std::string trajectoryID_v() const override {return "TimeAveragingMCWF_Trajectory";}

  const double relaxationTime_;
  
  mutable Averages averages_;
  mutable long sum_;
  
  const typename Averaged::Ptr av_;

};


} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_TIMEAVERAGINGMCWF_TRAJECTORY_H_INCLUDED
