// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_STOCHASTICTRAJECTORY_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_STOCHASTICTRAJECTORY_TCC_INCLUDED

#include "StochasticTrajectory.h"

#include "Conversions.h"
#include "Functional.h"
#include "ParsStochasticTrajectory.h"
#include "Range.h"
#include "impl/Trajectory.tcc"

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
  : Adaptive<A>(y,derivs,dtInit,p,scaleAbs,makerE),
    seed_(p.seed), isNoisy_(p.noise), randomized_(makerR(p.seed)) {}



template<typename T, typename T_ELEM>
std::ostream&
Ensemble<T,T_ELEM>::displayParameters_v(std::ostream& os) const
{
  return trajs_.begin()->displayParameters( os<<"# Ensemble of "<<trajs_.size()<<" trajectories."<<std::endl );
}


template<typename T, typename T_ELEM>
double
Ensemble<T,T_ELEM>::getDtDid_v() const
{
  return cpputils::accumulate(trajs_.begin(),trajs_.end(),0.,bind(&Trajectory::getDtDid,_1))/size2Double(trajs_.size());
}


template<typename T, typename T_ELEM>
void
Ensemble<T,T_ELEM>::evolve_v(double deltaT) const
{
  using namespace boost;

  if (log_) {
    progress_display pd(trajs_.size(),std::cerr);
    for (typename Impl::const_iterator i=trajs_.begin(); i!=trajs_.end(); ++i, ++pd) i->evolve(deltaT);
  }
  else
    for_each(trajs_,bind(&Trajectory::evolve,_1,deltaT));

}



template<typename T, typename T_ELEM>
const typename Ensemble<T,T_ELEM>::TBA_Type
Ensemble<T,T_ELEM>::averageInRange(size_t begin, size_t n) const
{
  return EnsembleTraits<T,T_ELEM>::averageInRange(trajs_.begin()+begin,trajs_.begin()+(begin+n),*this);

}




// Implementation of the traits class for the most commonly used case:

template<typename T, typename T_ELEM>
const typename EnsembleTraits<T,T_ELEM>::TBA_Type 
EnsembleTraits<T,T_ELEM>::averageInRange(typename Impl::const_iterator begin, typename Impl::const_iterator end, const ET&)
{
  return T(cpputils::accumulate(++begin,end,begin->toBeAveraged(),bind(&Elem::toBeAveraged,_1),cpputils::plus<T>())/size2Double(end-begin));
}


// A tentative specialization for the case when TBA_Type is a reference, it has a member function getInitializedTBA, and T_ELEM has a member function addTo

template<typename T, typename T_ELEM>
class EnsembleTraits<T&,T_ELEM>
{
public:
  typedef Ensemble<T&,T_ELEM> ET;

  typedef typename ET::Elem     Elem    ;
  typedef typename ET::Impl     Impl    ;
  typedef typename ET::TBA_Type TBA_Type;

  static const TBA_Type averageInRange(typename Impl::const_iterator begin, typename Impl::const_iterator end, const ET& et)
  {
    TBA_Type res(et.getInitializedTBA());
    
    for (typename Impl::const_iterator i=begin; i!=end; i++) i->toBeAveraged().addTo(res);

    return res/=size2Double(end-begin);

  }


};



} // trajectory


#endif // UTILS_INCLUDE_IMPL_STOCHASTICTRAJECTORY_TCC_INCLUDED
