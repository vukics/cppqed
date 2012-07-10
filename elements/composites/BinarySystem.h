// -*- C++ -*-
#ifndef _BINARY_SYSTEM_H
#define _BINARY_SYSTEM_H

#include "QuantumSystem.h"

#include "SubSystem.h"



namespace binary {

typedef structure::Interaction<2> Interaction;

typedef structure::SubSystemFree            SSF;
typedef structure::SubSystemsInteraction<2> SSI;



class Base
  : public structure::QuantumSystem<2>,
    public structure::Averaged     <2>
{
public:
  typedef structure::Averaged<1> Av1;
  typedef structure::Averaged<2> Av2;

  Base(const Interaction&);

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
  
};


#define CLASS_HEADER(Class) class Class : public structure::Class<2>

#define CLASS_BODY_PART(Class,Aux) public:				\
  typedef structure::Class<1> Aux##1;					\
  typedef structure::Class<2> Aux##2;					\
									\
  Class(const SSF& free0, const SSF& free1, const SSI& ia) : free0_(free0), free1_(free1), ia_(ia) {} \
									\
private:								\
  const SSF &free0_, &free1_;						\
  const SSI &ia_;							\



CLASS_HEADER(Exact)
{
  CLASS_BODY_PART(Exact,Ex)

  bool isUnitary() const;

  void actWithU(double, StateVectorLow&) const;

};


CLASS_HEADER(Hamiltonian)
{
  CLASS_BODY_PART(Hamiltonian,Ha)

  void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const;

};


CLASS_HEADER(Liouvillean)
{
  CLASS_BODY_PART(Liouvillean,Li)

  size_t              nJumps       ()                                   const;
  const Probabilities probabilities(double, const LazyDensityOperator&) const;
  void                actWithJ     (double, StateVectorLow&, size_t)    const;

  void displayKey(std::ostream&, size_t&) const;

};


#undef CLASS_BODY_PART
#undef CLASS_HEADER


} // binary



template<bool IS_EX=true, bool IS_HA=true, bool IS_LI=true>
class BinarySystem 
  : public binary::Base,
    public binary::Exact, 
    public binary::Hamiltonian,
    public binary::Liouvillean
{
public:
  typedef structure::Interaction<2> Interaction;

  BinarySystem(const Interaction&);

};



#endif // _BINARY_SYSTEM_H
