// -*- C++ -*-
#ifndef _STRUCTURE_LIOUVILLEAN_h
#define _STRUCTURE_LIOUVILLEAN_h

#include "LiouvilleanFwd.h"

#include "Types.h"

#include "Exception.h"


namespace structure {


class LiouvilleanCommon : public virtual LiouvilleanAveragedCommon
{
public:
  typedef TTD_DARRAY(1) Probabilities;
  // dpoverdt ladder type for Liouvilleans. Note that the extent is
  // not known at compile time because it depends on how many
  // subsystems there are which themselves are Liouvillean, which, in
  // turn, depends on parameters. This may also actually be only a
  // RANGE of a larger array if the present system is subsystem to a
  // larger system.

  static size_t nJumps(const LiouvilleanCommon* liouvillean)
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
class Liouvillean : public quantumdata::Types<RANK,LiouvilleanCommon>
{
public:
  typedef quantumdata::Types<RANK,LiouvilleanCommon> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;
  typedef typename Base::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename Base::Probabilities Probabilities;


  static const Probabilities probabilities(const LazyDensityOperator& matrix, const Liouvillean* liouvillean, StaticTag=theStaticOne) 
  {
    const Probabilities probas(liouvillean ? liouvillean->probabilities(matrix) : Probabilities());
#ifndef   NDEBUG
    if (size_t(probas.size())!=nJumps(liouvillean))
      throw LiouvilleanFishyException();
#endif // NDEBUG
    return probas;
  }


  static void actWithJ(StateVectorLow& psi, size_t jumpNo, const Liouvillean* liouvillean, StaticTag=theStaticOne)
  // jumpNo is the principal number of the jump to be performed
  {
    if (liouvillean) liouvillean->actWithJ(psi,jumpNo);
  }


  virtual ~Liouvillean() {}

  virtual void                actWithJ     (StateVectorLow&, size_t   ) const = 0;

private:    
  virtual const Probabilities probabilities(const LazyDensityOperator&) const = 0;

};


} // structure

#endif // _STRUCTURE_LIOUVILLEAN_h
