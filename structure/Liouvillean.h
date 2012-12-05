// -*- C++ -*-
#ifndef STRUCTURE_LIOUVILLEAN_H_INCLUDED
#define STRUCTURE_LIOUVILLEAN_H_INCLUDED

#include "LiouvilleanFwd.h"

#include "LiouvilleanAveragedCommon.h"


namespace structure {


template<int RANK>
class Liouvillean<RANK,true> : public quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> >
{
public:
  static const int N_RANK=RANK;

  typedef boost::shared_ptr<const Liouvillean> Ptr;

  typedef quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> > Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename Base::DArray1D Probabilities;

  virtual ~Liouvillean() {}

  void actWithJ(double t, StateVectorLow& psi, size_t jumpNo) const {return actWithJ_v(t,psi,jumpNo);}
  // jumpNo is the ordinal number of the jump to be performed

private:
  virtual void actWithJ_v(double, StateVectorLow&, size_t) const = 0;

};



template<int RANK>
class Liouvillean<RANK,false> : public Liouvillean<RANK,true>
{
public:
  typedef typename Liouvillean<RANK,true>::StateVectorLow      StateVectorLow     ;
  typedef typename Liouvillean<RANK,true>::LazyDensityOperator LazyDensityOperator;
  typedef typename Liouvillean<RANK,true>::Probabilities       Probabilities      ;

private:
  void                actWithJ_v(double, StateVectorLow& psi, size_t jumpNo) const {actWithJ_v(psi,jumpNo);}
  const Probabilities  average_v(double, const LazyDensityOperator&  matrix) const {return average_v(matrix);}

  virtual void                actWithJ_v(StateVectorLow&, size_t   ) const = 0;
  virtual const Probabilities  average_v(const LazyDensityOperator&) const = 0;

};



} // structure

#endif // STRUCTURE_LIOUVILLEAN_H_INCLUDED
