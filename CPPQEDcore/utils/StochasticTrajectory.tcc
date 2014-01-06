// -*- C++ -*-
#ifndef   UTILS_STOCHASTICTRAJECTORY_TCC_INCLUDED
#define   UTILS_STOCHASTICTRAJECTORY_TCC_INCLUDED

#include "StochasticTrajectory.h"

#include "Conversions.h"
#include "Functional.h"
#include "ParsStochasticTrajectory.h"
#include "Range.h"
#include "Trajectory.tcc"

#include <boost/progress.hpp>


namespace trajectory {


template<typename A, typename T> 
std::ostream&
Stochastic<A,T>::displayParameters_v(std::ostream& os) const
{
  return Adaptive<A>::displayParameters_v(os)<<"# Stochastic Trajectory Parameters: seed="<<seed_<<std::endl<<(isNoisy_ ? "" : "# No noise.\n");
}


template<typename A, typename T>
Stochastic<A,T>::Stochastic(A& y, typename Evolved::Derivs derivs,
                            double dtInit,
                            double epsRel, double epsAbs, const A& scaleAbs,
                            const evolved::Maker<A>& makerE,
                            unsigned long seed,
                            bool n,
                            const randomized::Maker& makerR)
  : Adaptive<A>(y,derivs,dtInit,epsRel,epsAbs,scaleAbs,makerE),
    seed_(seed), isNoisy_(n), randomized_(makerR(seed)) {}


template<typename A, typename T>
Stochastic<A,T>::Stochastic(A& y, typename Evolved::Derivs derivs,
                            double dtInit,
                            const A& scaleAbs,
                            const ParsStochastic& p,
                            const evolved::Maker<A>& makerE,
                            const randomized::Maker& makerR)
  : Stochastic(y,derivs,dtInit,p.epsRel,p.epsAbs,scaleAbs,makerE,p.seed,p.noise,makerR) {}



template<typename T, typename T_ELEM>
std::ostream&
Ensemble<T,T_ELEM>::displayParameters_v(std::ostream& os) const
{
  return trajs_.front().displayParameters( os<<"# Ensemble of "<<trajs_.size()<<" trajectories."<<std::endl );
}


template<typename T, typename T_ELEM>
double
Ensemble<T,T_ELEM>::getDtDid_v() const
{
  return cpputils::accumulate(trajs_.begin(),trajs_.end(),0.,bind(&Trajectory::getDtDid,_1))/size2Double(trajs_.size());
}


template<typename T, typename T_ELEM>
void
Ensemble<T,T_ELEM>::evolve_v(double deltaT)
{
  using namespace boost;

  if (log_) {
    progress_display pd(trajs_.size(),std::cerr);
    for (auto i=trajs_.begin(); i!=trajs_.end(); (++i, ++pd)) i->evolve(deltaT);
  }
  else
    for_each(trajs_,bind(&Trajectory::evolve,_1,deltaT));

}



template<typename T, typename T_ELEM>
auto
Ensemble<T,T_ELEM>::averageInRange(size_t begin, size_t n) const -> const ToBeAveragedType
{
  return ensemble::Traits<T,T_ELEM>::averageInRange(trajs_.begin()+begin,trajs_.begin()+(begin+n),*this);

}




// Implementation of the traits class for the most commonly used case:

template<typename T, typename T_ELEM>
auto
ensemble::Traits<T,T_ELEM>::averageInRange(typename Impl::const_iterator begin, typename Impl::const_iterator end, const EnsembleType&) -> const ToBeAveragedType
{
  return T(cpputils::accumulate(++begin,end,begin->toBeAveraged(),bind(&Elem::toBeAveraged,_1),cpputils::plus<T>())/size2Double(end-begin));
}


namespace ensemble {
// A tentative specialization for the case when ToBeAveragedType is a reference, it has a member function getInitializedToBeAveraged, and T_ELEM has a member function addTo

template<typename T, typename T_ELEM>
class Traits<T&,T_ELEM>
{
public:
  typedef Ensemble<T&,T_ELEM> EnsembleType;

  typedef typename EnsembleType::Elem             Elem            ;
  typedef typename EnsembleType::Impl             Impl            ;
  typedef typename EnsembleType::ToBeAveragedType ToBeAveragedType;

  static const ToBeAveragedType averageInRange(typename Impl::const_iterator begin, typename Impl::const_iterator end, const EnsembleType& et)
  {
    ToBeAveragedType res(et.getInitializedToBeAveraged());
    
    for (auto i=begin; i!=end; i++) i->toBeAveraged().addTo(res);

    return res/=size2Double(end-begin);

  }


};


} // ensemble


} // trajectory


#endif // UTILS_STOCHASTICTRAJECTORY_TCC_INCLUDED
