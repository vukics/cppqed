// -*- C++ -*-
#ifndef STRUCTURE_LIOUVILLEAN_H_INCLUDED
#define STRUCTURE_LIOUVILLEAN_H_INCLUDED

#include "LiouvilleanFwd.h"

#include "LiouvilleanAveragedCommon.h"
#include "Types.h"

#include "Exception.h"


namespace structure {


class LiouvilleanCommon : public LiouvilleanAveragedCommon
{
public:
  typedef boost::shared_ptr<const LiouvilleanCommon> Ptr;

  typedef DArray1D Probabilities;
  // dpoverdt ladder type for Liouvilleans. Note that the extent is not known at compile time because it depends on how many subsystems there are which themselves are Liouvillean, which, in turn, depends on parameters. This may also actually be only a RANGE of a larger array if the present system is subsystem to a larger system.

  virtual ~LiouvilleanCommon() {}

  size_t nJumps() const {return nJumps_v();}

private:
  virtual size_t nJumps_v() const = 0;

};


#ifndef   NDEBUG
struct LiouvilleanNumberMismatchException : cpputils::Exception {};
#endif // NDEBUG


template<int RANK>
class Liouvillean<RANK,true> : public quantumdata::Types<RANK,LiouvilleanCommon>
{
public:
  static const int N_RANK=RANK;

  typedef boost::shared_ptr<const Liouvillean> Ptr;

  typedef quantumdata::Types<RANK,LiouvilleanCommon> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename Base::Probabilities Probabilities;

  virtual ~Liouvillean() {}

  const Probabilities probabilities(double t, const LazyDensityOperator& matrix) const
  {
    const Probabilities probas(probabilities_v(t,matrix));
#ifndef   NDEBUG
    if (size_t(probas.size())!=Base::nJumps()) throw LiouvilleanNumberMismatchException();
#endif // NDEBUG
    return probas;
  }

  void actWithJ(double t, StateVectorLow& psi, size_t jumpNo) const {return actWithJ_v(t,psi,jumpNo);}
  // jumpNo is the ordinal number of the jump to be performed

private:
  virtual const Probabilities probabilities_v(double, const LazyDensityOperator&) const = 0;
  virtual void                     actWithJ_v(double, StateVectorLow&, size_t   ) const = 0;

};



template<int RANK>
class Liouvillean<RANK,false> : public Liouvillean<RANK,true>
{
public:
  typedef typename Liouvillean<RANK,true>::StateVectorLow      StateVectorLow     ;
  typedef typename Liouvillean<RANK,true>::LazyDensityOperator LazyDensityOperator;
  typedef typename Liouvillean<RANK,true>::Probabilities       Probabilities      ;

private:
  void                     actWithJ_v(double, StateVectorLow& psi, size_t jumpNo) const {actWithJ_v(psi,jumpNo);}
  const Probabilities probabilities_v(double, const LazyDensityOperator&  matrix) const {return probabilities_v(matrix);}

  virtual void                     actWithJ_v(StateVectorLow&, size_t   ) const = 0;
  virtual const Probabilities probabilities_v(const LazyDensityOperator&) const = 0;

};



} // structure

#endif // STRUCTURE_LIOUVILLEAN_H_INCLUDED
