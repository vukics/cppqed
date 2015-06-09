// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_TCC_INCLUDED

#include "StochasticTrajectory.h"

#include "Conversions.h"
#include "ParsStochasticTrajectory.h"
#include "Trajectory.tcc"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/adaptor/transformed.hpp>

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
  return accumulate(trajs_ | boost::adaptors::transformed(bind(&Trajectory::getDtDid,_1)) ,0.)/size2Double(trajs_.size());
}


template<typename T, typename T_ELEM>
void
Ensemble<T,T_ELEM>::evolve_v(double deltaT)
{
  using namespace boost;

  if (displayProgress_) {
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




// Naive generic implementation of the traits class
template<typename T, typename T_ELEM>
auto
ensemble::Traits<T,T_ELEM>::averageInRange(typename Impl::const_iterator begin, typename Impl::const_iterator end, const EnsembleType&) -> const ToBeAveragedType
{
  using namespace boost;
  return accumulate(make_iterator_range(++begin,end) | adaptors::transformed(bind(&Elem::toBeAveraged,_1)),begin->toBeAveraged())/size2Double(end-begin);
}


} // trajectory


#endif // CPPQEDCORE_UTILS_STOCHASTICTRAJECTORY_TCC_INCLUDED
