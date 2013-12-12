// -*- C++ -*-
#ifndef UTILS_EVOLVEDGSL_H_INCLUDED
#define UTILS_EVOLVEDGSL_H_INCLUDED

#include "Exception.h"

#include "Evolved.h"


namespace evolved {


class NonContiguousStorageException : public cpputils::Exception {};


template<typename A>
class MakerGSL : public Maker<A>
{
public:
  class _;
  
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


#endif // UTILS_EVOLVEDGSL_H_INCLUDED
