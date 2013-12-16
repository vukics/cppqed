/// \briefFile{Defines evolved::MakerGSL, which incorporates the actual \GSL-based implementation of evolved::Evolved}
// -*- C++ -*-
#ifndef UTILS_EVOLVEDGSL_H_INCLUDED
#define UTILS_EVOLVEDGSL_H_INCLUDED

#include "Exception.h"

#include "Evolved.h"


namespace evolved {


/// Thrown if the array supplied to MakerGSL has non contiguous storage
class NonContiguousStorageException : public cpputils::Exception {};


/// Implements Maker and incorporates MakerGSL::_, the actual \GSL-based implementation of the Evolved interface
template<typename A>
class MakerGSL : public Maker<A>
{
public:
  /**
   * \param nextDtTryCorrectionFactor This parameter is connected to the patching in _::step_v() of an undesired behaviour in `gsl_odeiv`: GSL seems to take dtDid as a basis of calculating dtTry,
   * instead of the previous dtTry. Now, since dtDid will be at most deltaT, this will result in a severe drop of dtTry, if deltaT is very small (e.g. when completing a given interval at its end).
   * Since this drop can be several orders of magnitude, this must not be allowed as it would slow down the simulation extremely. To patch this, we check before the actual timestep whether we are
   * in the case of a *very small* `deltaT` (`< dtTry/nextDtTryCorrectionFactor`). If this is found to be the case, then the present dtTry will be cached and used unchanged in the next step.
   * That is, the present (tiny) step is excluded from stepsize control.
   */
  MakerGSL(SteppingFunction sf=SF_RKCK, double nextDtTryCorrectionFactor=100.) : sf_(sf), nextDtTryCorrectionFactor_(nextDtTryCorrectionFactor) {}

private:
  class _;
  
  typedef typename Maker<A>::Ptr Ptr;
  typedef typename Maker<A>::Derivs Derivs;

  const Ptr make(A&, Derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs) const;

  const SteppingFunction sf_;
  const double nextDtTryCorrectionFactor_;
  
};


namespace details {

class Impl;

typedef boost::shared_ptr<Impl> ImplPtr;

extern const int onSuccess;

} // details

} // evolved


#endif // UTILS_EVOLVEDGSL_H_INCLUDED
