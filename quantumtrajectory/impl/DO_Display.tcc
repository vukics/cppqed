// -*- C++ -*-
#ifndef   _DO_DISPLAY_IMPL_H
#define   _DO_DISPLAY_IMPL_H

#include "DimensionsBookkeeper.h"
#include "NegPT.h"
#include "Structure.h"

#include "FormDouble.h"
#include "ParsTrajectory.h"

#include<string>
#include<sstream>


namespace quantumtrajectory {


namespace details {


template<int RANK, typename V>
DO_Display<RANK,V>::DO_Display(const QuantumSystem& qs,
				     const ParsTrajectory& p,
				     bool negativity,
				     size_t equalCount) throw(DimensionalityMismatchException)
  : av_(structure::qsa(&qs)),
    negativity_(negativity),    
    autoStop_(p.autoStop),
    lastCrit_(),
    equalCount_(equalCount)
{
}



template<int RANK, typename V>
size_t
DO_Display<RANK,V>::displayMoreKey(std::ostream& os) const 
{
  size_t i=3; 
  Averaged::displayKey(os,i,av_); 
  if (negativity_) os<<"# Trajectory "<<i<<". negativity"<<std::endl;
  return i;
}


template<int RANK, typename V>
void
DO_Display<RANK,V>::displayMore(double t, const DensityOperator& rho, std::ostream& os, int precision) const 
  throw(StoppingCriterionReachedException)
{
  using namespace std;
  stringstream line(stringstream::in | stringstream::out);
  {
    if (av_) 
      Averaged::display(t,rho,line,precision,av_);
    if (negativity_) line<<'\t'<<formdouble::FormDouble(precision)(quantumdata::negPT(rho,V()));
    line<<endl;
  }

  // A heuristic stopping criterion: If in column autoStop the same
  // value (in the given precision) is Displayed equalCountLimit_
  // times in a row, then stop.
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

}


} // details


} // quantumtrajectory


#endif // _DO_DISPLAY_IMPL_H
