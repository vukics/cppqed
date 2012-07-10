// -*- C++ -*-
#ifndef _BINARY_SYSTEM_H
#define _BINARY_SYSTEM_H

#include "QuantumSystem.h"

#include "SubSystem.h"


namespace binary {

typedef blitz::TinyVector<bool,3> Mask;

typedef structure::Interaction<2> Interaction;

typedef structure::SubSystemFree            SSF;
typedef structure::SubSystemsInteraction<2> SSI;


class Base
  : public structure::QuantumSystem<2>,
    public structure::Averaged<2>
{
public:
  typedef structure::Averaged<1> Av1;
  typedef structure::Averaged<2> Av2;

  Base(const Interaction&, const Mask&);

  const SSF& getFree0() const {return free0_;}
  const SSF& getFree1() const {return free1_;}
  const SSI& getIA   () const {return    ia_;}

private:
  double highestFrequency (             ) const;
  void   displayParameters(std::ostream&) const;

  void   displayKey(std::ostream&, size_t&) const;
  size_t nAvr      (                      ) const;

  const Averages average(double, const LazyDensityOperator&)  const;
  void           process(Averages&)                           const;
  void           display(const Averages&, std::ostream&, int) const;

  const SSF free0_, free1_;

  const SSI ia_;
  
  const Mask mask_;

};



class Liouvillean : public structure::Liouvillean<2>
{
public:
  typedef structure::Liouvillean  <1> Li1; 
  typedef structure::Liouvillean  <2> Li2; 

  Liouvillean(const SSF&, const SSF&, const SSI&, const Mask&);

private:
  size_t              nJumps       ()                                   const;
  const Probabilities probabilities(double, const LazyDensityOperator&) const;
  void                actWithJ     (double, StateVectorLow&, size_t)    const;

  void displayKey(std::ostream&, size_t&) const;

  const SSF &free0_, &free1_;
  const SSI &ia_;

  const Mask mask_;

};



} // binary



class BinarySystem 
  : public binary::Base,
    public structure::Exact        <2>, 
    public structure::Hamiltonian  <2>,
    public binary::Liouvillean
{
public:
  typedef structure::Interaction<2> Interaction;
  typedef quantumdata::Types<2> Types;

  typedef Types::    StateVectorLow     StateVectorLow;
  typedef Types::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<2> LazyDensityOperator;

  typedef structure::Exact        <1> Ex1; 
  typedef structure::Hamiltonian  <1> Ha1;
  typedef structure::Liouvillean  <1> Li1; 

  typedef structure::Exact        <2> Ex2; 
  typedef structure::Hamiltonian  <2> Ha2;
  typedef structure::Liouvillean  <2> Li2; 

  BinarySystem(const Interaction&);

private:
  bool isUnitary() const;

  void actWithU(double, StateVectorLow&) const;

  void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const;

  const binary::SSF free0_, free1_;
  const binary::SSI ia_;

};



#endif // _BINARY_SYSTEM_H
