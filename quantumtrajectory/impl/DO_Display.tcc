// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_DO_DISPLAY_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_DO_DISPLAY_TCC_INCLUDED

#include "DO_Display.h"

#include "DimensionsBookkeeper.h"
#include "impl/NegPT.tcc"
#include "Structure.h"

#include "FormDouble.h"

#include <string>
#include <sstream>


namespace quantumtrajectory {


namespace details {


template<int RANK, typename V>
DO_Display<RANK,V>::DO_Display(AveragedPtr av,
                               const ParsEvolved& p,
                               bool negativity,
                               size_t equalCount) throw(DimensionalityMismatchException)
  : av_(av),
    negativity_(negativity),    
    autoStop_(0), // Disabled at the moment
    lastCrit_(),
    equalCount_(equalCount)
{
}



template<int RANK, typename V>
std::ostream&
DO_Display<RANK,V>::displayKey(std::ostream& os, size_t& i) const 
{
  if (av_) av_->displayKey(os,i); 
  if (negativity_) os<<"# Trajectory "<<i<<". negativity"<<std::endl;
  return os;
}


template<int RANK, typename V>
std::ostream&
DO_Display<RANK,V>::display(double t, const DensityOperator& rho, std::ostream& os, int precision) const 
  throw(StoppingCriterionReachedException)
{
  using namespace std;
  stringstream line(stringstream::in | stringstream::out);
  {
    structure::display(av_,t,rho,line,precision);
    if (negativity_) line<<'\t'<<FormDouble(precision)(quantumdata::negPT(rho,V()));
  }

  // A heuristic stopping criterion: If in column autoStop the same value (in the given precision) is Displayed equalCountLimit_ times in a row, then stop.
  if (autoStop_ && av_) {
    double crit;
    for (size_t i=0; i<autoStop_; i++)
      line>>crit;
    if (crit==lastCrit_) --equalCount_;
    else equalCount_=0;
    lastCrit_=crit;
  }

  os<<line.str();
  if (!equalCount_) throw StoppingCriterionReachedException();

  return os;
}


} // details


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_IMPL_DO_DISPLAY_TCC_INCLUDED
