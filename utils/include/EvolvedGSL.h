// -*- C++ -*-
#ifndef _EVOLVED_GSL_H
#define _EVOLVED_GSL_H

#include "EvolvedGSLFwd.h"

#include "Exception.h"

#include "Evolved.h"


namespace evolved {


template<typename A>
class MakerGSL : public Maker<A>
{
public:
  typedef typename Maker<A>::SmartPtr SmartPtr;
  typedef typename Maker<A>::Derivs   Derivs  ;

  const SmartPtr operator()(A&, Derivs, double dtInit, double epsRel, double epsAbs, const A& scaleAbs) const;
  
};



namespace details  {


///////////////////////
//
// Forward declarations 
//
///////////////////////

typedef boost::shared_ptr<Impl> ImplSmartPtr;

// Two indirections needed because cannot declare member functions for Impl here

ImplSmartPtr createImpl(void*, size_t, int(double,const double*,double*,void*), double, double, const double*);

void apply(ImplSmartPtr, double*, double, double*, double*);

extern const int onSuccess;

class NonContiguousStorageException : public cpputils::Exception {};


/////////////////////////////
//
// Evolved GSL implementation
//
/////////////////////////////

template<typename A>
class EvolvedGSL : public Evolved<A> 
{
public:
  typedef Evolved<A> Base;

  typedef typename Base::Derivs Derivs;
  typedef typename Base::Traits Traits;

  using Base::getA; using Base::getTime; using Base::getDtTry;

  EvolvedGSL(
	     A&,
	     Derivs,
	     double,
	     double,
	     double,
	     const A&
	     );

  void step(double);
  
private:
  ImplSmartPtr pImpl_;

};


} // details


} // evolved


#include "impl/EvolvedGSL.tcc"


#endif // _EVOLVED_GSL_H
