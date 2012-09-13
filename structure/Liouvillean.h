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

  static size_t nJumps(Ptr liouvillean)
  {
    return liouvillean ? liouvillean->nJumps() : 0;
  }

  virtual ~LiouvilleanCommon() {}

private:
  virtual size_t nJumps() const = 0;

};


#ifndef   NDEBUG
struct LiouvilleanFishyException : cpputils::Exception {};
#endif // NDEBUG


template<int RANK>
class Liouvillean<RANK,true> : public quantumdata::Types<RANK,LiouvilleanCommon>
{
public:
  typedef boost::shared_ptr<const Liouvillean> Ptr;

  typedef quantumdata::Types<RANK,LiouvilleanCommon> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename Base::Probabilities Probabilities;


  static const Probabilities probabilities(double t, const LazyDensityOperator& matrix, Ptr liouvillean, StaticTag=theStaticOne)
  {
    const Probabilities probas(liouvillean ? liouvillean->probabilities(t,matrix) : Probabilities());
#ifndef   NDEBUG
    if (size_t(probas.size())!=Base::nJumps(liouvillean))
      throw LiouvilleanFishyException();
#endif // NDEBUG
    return probas;
  }


  static void actWithJ(double t, StateVectorLow& psi, size_t jumpNo, Ptr liouvillean, StaticTag=theStaticOne)
  // jumpNo is the ordinal number of the jump to be performed
  {
    if (liouvillean) liouvillean->actWithJ(t,psi,jumpNo);
  }


  virtual ~Liouvillean() {}

  virtual void                actWithJ     (double, StateVectorLow&, size_t   ) const = 0;

private:    
  virtual const Probabilities probabilities(double, const LazyDensityOperator&) const = 0;

};



template<int RANK>
class Liouvillean<RANK,false> : public Liouvillean<RANK,true>
{
public:
  typedef typename Liouvillean<RANK,true>::StateVectorLow      StateVectorLow     ;
  typedef typename Liouvillean<RANK,true>::LazyDensityOperator LazyDensityOperator;
  typedef typename Liouvillean<RANK,true>::Probabilities       Probabilities      ;


private:
  void                actWithJ     (double, StateVectorLow& psi, size_t jumpNo) const {actWithJ(psi,jumpNo);}
  const Probabilities probabilities(double, const LazyDensityOperator&  matrix) const {return probabilities(matrix);}

  virtual void                actWithJ     (StateVectorLow&, size_t   ) const = 0;
  virtual const Probabilities probabilities(const LazyDensityOperator&) const = 0;

};



} // structure

#endif // STRUCTURE_LIOUVILLEAN_H_INCLUDED
