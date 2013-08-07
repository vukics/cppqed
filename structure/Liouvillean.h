// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_LIOUVILLEAN_H_INCLUDED
#define STRUCTURE_LIOUVILLEAN_H_INCLUDED

#include "LiouvilleanFwd.h"

#include "LiouvilleanAveragedCommon.h"


namespace structure {


/// The first partial specialization of the general template Liouvillean for the one-time dependence case (\link TimeDependence Cases 1 & 2\endlink)
template<int RANK>
class Liouvillean<RANK,true> : public quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> >
{
public:
  static const int N_RANK=RANK;

  typedef boost::shared_ptr<const Liouvillean> Ptr;

private:
  typedef quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> > Base;

public:
  typedef typename Base::    StateVectorLow     StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::DArray1D Probabilities; ///< The 1D real array for storing the jump probabilities

  virtual ~Liouvillean() {}

  /// Performs the quantum jump operation \f$\ket\Psi\rightarrow J_m(t)\ket\Psi\f$
  void actWithJ(double t,            ///<[in] \f$t\f$
                StateVectorLow& psi, ///<[in/out] \f$\ket\Psi\f$
                size_t m             ///<[out] \f$m\f$
               ) const {return actWithJ_v(t,psi,m);}

private:
  virtual void actWithJ_v(double, StateVectorLow&, size_t) const = 0;

};


/// The first partial specialization of the general template Liouvillean for the no-time dependence case (\link TimeDependence Cases 3 & 4\endlink)
template<int RANK>
class Liouvillean<RANK,false> : public Liouvillean<RANK,true>
{
public:
  typedef typename Liouvillean<RANK,true>::StateVectorLow      StateVectorLow     ;
  typedef typename Liouvillean<RANK,true>::LazyDensityOperator LazyDensityOperator;
  typedef typename Liouvillean<RANK,true>::Probabilities       Probabilities      ;

private:
  void                actWithJ_v(double, StateVectorLow& psi, size_t jumpNo) const {actWithJ_v(psi,jumpNo);}   ///< Redirects the virtual inherited from Liouvillean<RANK,true>
  const Probabilities  average_v(double, const LazyDensityOperator&  matrix) const {return average_v(matrix);} ///< Redirects the virtual inherited from LiouvilleanAveragedCommonRanked

  virtual void                actWithJ_v(StateVectorLow&, size_t   ) const = 0;
  virtual const Probabilities  average_v(const LazyDensityOperator&) const = 0;

};



} // structure

#endif // STRUCTURE_LIOUVILLEAN_H_INCLUDED
