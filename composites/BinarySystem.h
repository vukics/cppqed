// -*- C++ -*-
#ifndef _BINARY_SYSTEM_H
#define _BINARY_SYSTEM_H

#include "QuantumSystem.h"

#include "SubSystem.h"

typedef structure::Interaction<2> Interaction;


class BinarySystem 
  : public structure::QuantumSystem<2>, 
    public structure::Exact        <2>, 
    public structure::Hamiltonian  <2>,
    public structure::Liouvillean  <2>,
    public structure::Averaged     <2>
{
public:
  typedef quantumdata::Types<2> Types;

  typedef Types::    StateVectorLow     StateVectorLow;
  typedef Types::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<2> LazyDensityOperator;

  typedef structure::QuantumSystem<1> QS1;
  typedef structure::Exact        <1> Ex1; 
  typedef structure::Hamiltonian  <1> Ha1;
  typedef structure::Liouvillean  <1> Li1; 
  typedef structure::Averaged     <1> Av1;

  typedef structure::QuantumSystem<2> QS2;
  typedef structure::Exact        <2> Ex2; 
  typedef structure::Hamiltonian  <2> Ha2;
  typedef structure::Liouvillean  <2> Li2; 
  typedef structure::Averaged     <2> Av2;

  BinarySystem(const Interaction&);

private:
  double highestFrequency (             ) const;
  void   displayParameters(std::ostream&) const;
  
  bool isUnitary() const;

  void actWithU(double, StateVectorLow&) const;

  void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const;

  size_t              nJumps       ()                           const;
  const Probabilities probabilities(const LazyDensityOperator&) const;
  void                actWithJ     (StateVectorLow&, size_t)    const;

  void   displayKey(std::ostream&, size_t&) const;
  size_t nAvr      ()                       const;

  const Averages average(const LazyDensityOperator&)          const;
  void           process(Averages&)                           const;
  void           display(const Averages&, std::ostream&, int) const;


  const structure::SubSystemFree free0_;
  const structure::SubSystemFree free1_;

  const structure::SubSystemsInteraction<2> ia_;

};



#endif // _BINARY_SYSTEM_H
