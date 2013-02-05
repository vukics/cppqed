// -*- C++ -*-
#ifndef UTILS_INCLUDE_EVOLVEDGSL_H_INCLUDED
#define UTILS_INCLUDE_EVOLVEDGSL_H_INCLUDED

#include "Exception.h"

#include "Evolved.h"


namespace evolved {


enum SteppingFunction {SF_RKCK, SF_RK8PD};

std::ostream& operator<<(std::ostream&, SteppingFunction);
std::istream& operator>>(std::istream&, SteppingFunction&);


class NonContiguousStorageException : public cpputils::Exception {};


template<typename A>
class MakerGSL : public Maker<A>
{
public:
  class GSL;
  
  typedef typename Maker<A>::Ptr Ptr;
  typedef typename Maker<A>::Derivs   Derivs  ;

  MakerGSL(SteppingFunction sf=SF_RKCK, double nextDtTryCorretionFactor=100.) : sf_(sf), nextDtTryCorretionFactor_(nextDtTryCorretionFactor) {}

  const Ptr operator()(A&, Derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs) const;

private:
  const SteppingFunction sf_;
  const double nextDtTryCorretionFactor_;
  
};


namespace details {

class Impl;

typedef boost::shared_ptr<Impl> ImplPtr;

extern const int onSuccess;

} // details

} // evolved


#endif // UTILS_INCLUDE_EVOLVEDGSL_H_INCLUDED
