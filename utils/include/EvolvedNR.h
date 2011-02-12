// -*- C++ -*-
#ifndef _EVOLVED_NR_H
#define _EVOLVED_NR_H

#include "EvolvedNRFwd.h"

#include "Evolved.h"


namespace evolved {


template<typename A> 
class MakerNR : public Maker<A> 
{
public:
  typedef typename Maker<A>::SmartPtr SmartPtr;
  typedef typename Maker<A>::Derivs   Derivs  ;

  const SmartPtr operator()(A&, Derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs) const;
  
};



namespace details {

// Note that EvolvedNR assumes that Array features vector-space
// operations. This could be released by the use of overloaded worker
// functions.

///////////////////////////////////
//
// Evolved Num. Rec. implementation
//
///////////////////////////////////

template<typename A> 
class EvolvedNR : public Evolved<A>
{
public:
  typedef Evolved<A> Base;

  typedef typename Base::Derivs Derivs;
  typedef typename Base::Traits Traits;

  using Base::getA; using Base::getTime; using Base::getDtDid; using Base::getDtTry;

  EvolvedNR(A&, Derivs, double, double, double, const A&);

  void step(double);

private:
  A yscal_, dydt_;
  const A& scaleAbs_;

};


} // details


} // evolved

#include "impl/EvolvedNR.tcc"


#endif // _EVOLVED_NR_H
