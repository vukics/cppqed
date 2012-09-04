// -*- C++ -*-
#ifndef STRUCTURE_INTERACTION_H_INCLUDED
#define STRUCTURE_INTERACTION_H_INCLUDED

/*
  The simplest possible implementation allowing for interaction
  between frees only. If we are to achieve recursiveness, that is the
  possibility of nesting composite systems into even more composite
  ones, we have to allow interaction between composites as well. So
  QuantumSystems in general.

  For this, Interaction should be even more templated taking
  compile-time vectors. These specify between which quantum numbers of
  the subsystems the interaction acts. Obviously, as many compile-time
  vectors are needed as the number of subsystems.
*/

#include "InteractionFwd.h"

#include "DynamicsBase.h"
#include "Free.h"

#include "BlitzTiny.h"


namespace structure {


template<int RANK>
class Interaction : public DynamicsBase
{
public:
  typedef blitz::TinyVector<Free::SmartPtr,RANK> Frees;

  explicit Interaction(const Frees& frees, 
		       const    RealFreqs&    realFreqs=   RealFreqs(), 
		       const ComplexFreqs& complexFreqs=ComplexFreqs())
    : DynamicsBase(realFreqs,complexFreqs), frees_(frees) {}

  const Frees& getFrees() const {return frees_;}

private:
  const Frees frees_;

};


} // structure

#endif // STRUCTURE_INTERACTION_H_INCLUDED
