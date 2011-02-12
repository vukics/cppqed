// -*- C++ -*-
#ifndef   _STOCHASTIC_TRAJECTORY_IMPL_H
#define   _STOCHASTIC_TRAJECTORY_IMPL_H

#include "Algorithm.h"
#include "Functional.h"

#include "Range.h"

#include<boost/bind.hpp>


namespace trajectory {


template<typename A, typename T> 
void 
StochasticTrajectory<A,T>::displayParameters() const
{
  Trajectory<A>::displayParameters();
  TrajectoryBase::getOstream()<<"# Stochastic Trajectory Parameters: seed="<<seed_<<std::endl;
  if (!isNoisy_) 
    TrajectoryBase::getOstream()<<"# No noise."<<std::endl;
}


template<typename A, typename T>
StochasticTrajectory<A,T>::StochasticTrajectory(A& y, typename Evolved::Derivs derivs,
						double dtInit,
						double epsRel, double epsAbs, const A& scaleAbs,
						const evolved::Maker<A>& makerE,
						unsigned long seed,
						bool n,
						const randomized::Maker& makerR)
  : Trajectory<A>(y,derivs,dtInit,epsRel,epsAbs,scaleAbs,makerE),
    seed_(seed), isNoisy_(n), randomized_(makerR(seed)) {}


template<typename A, typename T>
StochasticTrajectory<A,T>::StochasticTrajectory(A& y, typename Evolved::Derivs derivs,
						double dtInit,
						const A& scaleAbs,
						const ParsStochasticTrajectory& p,
						const evolved::Maker<A>& makerE,
						const randomized::Maker& makerR)
  : Trajectory<A>(y,derivs,dtInit,scaleAbs,p,makerE),
    seed_(p.seed), isNoisy_(p.noise), randomized_(makerR(p.seed)) {}



template<typename T, typename T_ELEM>
void 
EnsembleTrajectories<T,T_ELEM>::displayParameters() const
{
  TrajectoryBase::getOstream()<<"# Ensemble of "<<trajs_.size()<<" trajectories."<<std::endl;
  trajs_.begin()->displayParameters();
}


template<typename T, typename T_ELEM>
double
EnsembleTrajectories<T,T_ELEM>::getDtDid() const
{
  return cpputils::accumulate(trajs_.begin(),trajs_.end(),0.,bind(&TrajectoryBase::getDtDid,_1))/size2Double(trajs_.size());
}


template<typename T, typename T_ELEM>
void
EnsembleTrajectories<T,T_ELEM>::evolve(double deltaT) const
{
  boost::for_each(trajs_,bind(&TrajectoryBase::evolve,_1,deltaT));
}



template<typename T, typename T_ELEM>
const typename EnsembleTrajectories<T,T_ELEM>::TBA_Type
EnsembleTrajectories<T,T_ELEM>::toBeAveraged() const
{
  return EnsembleTrajectoriesTraits<T,T_ELEM>::toBeAveraged(*this);

}




// Implementation of the traits class for the most commonly used case:

template<typename T, typename T_ELEM>
const typename EnsembleTrajectoriesTraits<T,T_ELEM>::TBA_Type 
EnsembleTrajectoriesTraits<T,T_ELEM>::toBeAveraged(const ET& et)
{
  const Impl& trajs(et.getTrajs());
  return T(cpputils::accumulate(++trajs.begin(),trajs.end(),trajs.begin()->toBeAveraged(),
				bind(&Elem::toBeAveraged,_1),
				cpputils::plus<T>())/size2Double(trajs.size()));
}


// A tentative specialization for the case when IS_REF=true and T_ELEM has a member function addTo

template<typename T, typename T_ELEM>
class EnsembleTrajectoriesTraits<T&,T_ELEM>
{
public:
  typedef EnsembleTrajectories<T&,T_ELEM> ET;

  typedef typename ET::Elem     Elem    ;
  typedef typename ET::Impl     Impl    ;
  typedef typename ET::TBA_Type TBA_Type;

  static const TBA_Type toBeAveraged(const ET& et)
  {
    TBA_Type res(et.getInitializedTBA());
    const Impl& trajs(et.getTrajs());
    
    for (typename Impl::const_iterator i=trajs.begin(); i!=trajs.end(); i++) i->toBeAveraged().addTo(res);

    return res/=size2Double(trajs.size());

  }


};



} // trajectory


#endif // _STOCHASTIC_TRAJECTORY_IMPL_H
