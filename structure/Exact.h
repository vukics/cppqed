// -*- C++ -*-
#ifndef STRUCTURE_EXACT_H_INCLUDED
#define STRUCTURE_EXACT_H_INCLUDED

#include "ExactFwd.h"

#include "Types.h"


namespace structure {


class ExactCommon
{
public:
  static  bool isUnitary(const ExactCommon* exactCommon) {return exactCommon ? exactCommon->isUnitary() : true;}

  virtual ~ExactCommon() {}

private:
  virtual bool isUnitary() const = 0;

};


template<int RANK>
class Exact : public ExactCommon, private quantumdata::Types<RANK> 
{
public:
  typedef quantumdata::Types<RANK> Base;
  typedef typename Base::StateVectorLow StateVectorLow;

  static  void actWithU(double t, StateVectorLow& psi, const Exact* exact, StaticTag=theStaticOne) {if (exact) exact->actWithU(t,psi);}
  // The exact (in general, non-unitary) part of evolution. Put
  // otherwise, the operator transfoming between normal and
  // interaction pictures.

  // The tag is for helping boost.bind to distinguish between this
  // function and the virtual one below. Simply permuting the
  // arguments is not enough for this!

  virtual ~Exact() {}

  virtual void actWithU(double, StateVectorLow&) const = 0;

};


} // structure

#endif // STRUCTURE_EXACT_H_INCLUDED
