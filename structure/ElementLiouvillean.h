// -*- C++ -*-
#ifndef _ELEMENT_H
#define _ELEMENT_H

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"

#ifndef   NDEBUG
#include "Exception.h"
#endif // NDEBUG

#include<boost/function.hpp>


namespace structure {



template<int RANK, int NOJ>
class ElementLiouvillean : public Liouvillean<RANK>
{
public:
  typedef Liouvillean<RANK> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename LiouvilleanCommon::Probabilities Probabilities;

  typedef boost::function<void  (      StateVectorLow&     )> JumpStrategy;
  typedef boost::function<double(const LazyDensityOperator&)> JumpProbabilityStrategy;

  typedef blitz::TinyVector<JumpStrategy           ,NOJ> JumpStrategies;
  typedef blitz::TinyVector<JumpProbabilityStrategy,NOJ> JumpProbabilityStrategies;

protected:
  ElementLiouvillean(const JumpStrategies& jumps, const JumpProbabilityStrategies& jumpProbas) 
    : jumps_(jumps), jumpProbas_(jumpProbas) {}

private:
  size_t nJumps() const {return NOJ;}

  const Probabilities probabilities(const LazyDensityOperator&) const;

  void actWithJ(StateVectorLow& psi, size_t jumpNo) const {jumps_(jumpNo)(psi);}

  const JumpStrategies            jumps_     ;
  const JumpProbabilityStrategies jumpProbas_;

};


template<int RANK>
class ElementLiouvillean<RANK,1> : public Liouvillean<RANK>
// This specialization can use the virtual-function technique of old
{
public:
  typedef Liouvillean<RANK> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename LiouvilleanCommon::Probabilities Probabilities;

private:
  size_t nJumps() const {return 1;}

  const Probabilities probabilities(const LazyDensityOperator&) const;

#ifndef   NDEBUG
  struct ElementLiouvilleanException : cpputils::Exception {};
#endif // NDEBUG

  void actWithJ(StateVectorLow& psi, size_t 
#ifndef   NDEBUG
		  jumpNo
#endif // NDEBUG
		  ) const {
#ifndef   NDEBUG
    if (jumpNo) throw ElementLiouvilleanException(); 
#endif // NDEBUG
    doActWithJ(psi);
  }


  virtual void   doActWithJ (      StateVectorLow     &) const = 0;
  virtual double probability(const LazyDensityOperator&) const = 0;

};


} // structure

#include "impl/ElementLiouvillean.tcc"

#endif // _ELEMENT_H
