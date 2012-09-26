// -*- C++ -*-
#ifndef STRUCTURE_EXACT_H_INCLUDED
#define STRUCTURE_EXACT_H_INCLUDED

#include "ExactFwd.h"

#include "Types.h"

#include <boost/shared_ptr.hpp>

namespace structure {


class ExactCommon
{
public:
  typedef boost::shared_ptr<const ExactCommon> Ptr;

  virtual ~ExactCommon() {}

  bool isUnitary() const {return isUnitary();}

private:
  virtual bool isUnitary_v() const = 0;

};


template<int RANK>
class Exact : public ExactCommon, private quantumdata::Types<RANK> 
{
public:
  typedef boost::shared_ptr<const Exact> Ptr;

  typedef quantumdata::Types<RANK> Base;
  typedef typename Base::StateVectorLow StateVectorLow;

  static  void actWithU(double t, StateVectorLow& psi, Ptr exact, StaticTag=theStaticOne) {if (exact) exact->actWithU(t,psi);}
  // The exact (in general, non-unitary) part of evolution. Put otherwise, the operator transfoming between normal and interaction pictures.

  virtual ~Exact() {}

  virtual void actWithU(double, StateVectorLow&) const = 0;

};


} // structure

#endif // STRUCTURE_EXACT_H_INCLUDED
