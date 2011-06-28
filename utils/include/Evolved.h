// -*- C++ -*-
/*
  Evolved is intended as a common interface for (adaptive stepsize)
  ODE drivers. It takes the Array it operates on as template
  parameter. A given Array can be adapted to the form expected by
  Evolved by a suitable specialization of ArrayMemoryTraits. 

  The Array which is actually "Evolved" is by no means owned by
  Evolved (meaning that it is not deallocated when an Evolved is
  destructed).

  Note that the function calculating the derivative cannot be
  processed as metadata: Indeed, usually it's a functor, which
  incorporates additional data obtained at run-time.
*/

#ifndef _EVOLVED_H
#define _EVOLVED_H

#include "EvolvedFwd.h"

#include "ArrayTraitsFwd.h"

#include<boost/shared_ptr.hpp> // instead of std::tr1::shared_ptr
#include<boost/function.hpp>   // instead of std::tr1::function
#include<boost/utility.hpp>

namespace evolved {

  
namespace details {


class EvolvedCommon : private boost::noncopyable 
{
public:
  double getTime(        ) const {return t_;}
  void   setTime(double t)       {t_=t;}

  double getDtDid() const {return dtDid_;}

  double getDtTry(            ) const {return dtTry_;}
  void   setDtTry(double dtTry)       {dtTry_=dtTry;}

  double getEpsRel() const {return epsRel_;}
  double getEpsAbs() const {return epsAbs_;}

  void update(double t, double dtTry) {dtDid_=t-t_; t_=t; dtTry_=dtTry;}

protected:
  EvolvedCommon(double dtInit, double epsRel, double epsAbs);

  void setDtDid(double dtDid) {dtDid_=dtDid;}

private:
  double t_, dtTry_, dtDid_;

  const double epsRel_, epsAbs_;

};

} // details


////////////////////
//
// Evolved interface
//
////////////////////

template<typename A>
class Evolved : public details::EvolvedCommon
{
public:
  typedef cpputils::ArrayMemoryTraits<A> Traits;

  typedef boost::function<void(double, const A&, A&)> Derivs;

  typedef boost::shared_ptr<Evolved> SmartPtr;


  Evolved(A&, Derivs, double dtInit, double epsRel, double epsAbs);

  virtual ~Evolved() {}
    
  virtual void step(double deltaT) = 0;
  // Takes a single adaptive step of maximum length deltaT

        A& getA()       {return a_;}
  const A& getA() const {return a_;}

  const Derivs getDerivs() const {return derivs_;}

private:
  A& a_;

  Derivs derivs_;

};


template<typename E>
void evolve(E&, double deltaT);
// evolves for exactly deltaT


////////////////
//
// factory class
//
////////////////

template<typename A>
class Maker
{
public:
  typedef typename Evolved<A>::SmartPtr SmartPtr;
  typedef typename Evolved<A>::Derivs   Derivs  ;
  
  virtual const SmartPtr operator()(A&, Derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs) const = 0;

  virtual ~Maker() {}
  
};



} // evolved

#include "impl/Evolved.tcc"

/*
The problem with the original solution of defining a static factory
function inside Evolved together with an enum for the implementations
is that all the implementation templates need to be instantiated,
which may not be possible for all Array types.

E.g.

template<typename A>
typename Evolved<A>::SmartPtr Evolved<A>::create(
						 A& b,
						 Derivs derivs,
						 double dtInit,
						 double epsRel,
						 double epsAbs,
						 const A& scaleAbs,
						 Implementation ei
						 )
{
  if (ei==NR) 
    return typename Evolved<A>::SmartPtr(new details:: NR_Evolved<A>(b,Derivs,dtInit,epsRel,epsAbs,scaleAbs));
  else
    return typename Evolved<A>::SmartPtr(new details::GSL_Evolved<A>(b,Derivs,dtInit,epsRel,epsAbs,scaleAbs));
}


The factory CLASS technique, however, allows for both compile-time &
run-time implementation-selection, and also for implementations whose
constructors take additional arguments.

*/



#endif // _EVOLVED_H
