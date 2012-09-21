// -*- C++ -*-
#ifndef ELEMENTS_COMPOSITES_BINARYSYSTEM_H_INCLUDED
#define ELEMENTS_COMPOSITES_BINARYSYSTEM_H_INCLUDED

#include "BinarySystemFwd.h"

#include "QuantumSystem.h"
#include "SubSystem.h"

#include "SmartPtr.h"


namespace binary {


typedef boost::shared_ptr<const Base> SmartPtr;

typedef structure::Interaction<2> Interaction;


const SmartPtr make(Interaction::Ptr);

template<typename IA>
const SmartPtr make(const IA& ia)
{
  return make(cpputils::sharedPointerize(ia));
}


typedef composite::SubSystemFree            SSF;
typedef composite::SubSystemsInteraction<2> SSI;



class Base
  : public structure::QuantumSystem<2>,
    public structure::Averaged     <2>
{
public:
  typedef structure::Averaged<1> Av1;
  typedef structure::Averaged<2> Av2;

  Base(Interaction::Ptr);

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


template<typename>
class EmptyBase
{
public:
  EmptyBase(const SSF&, const SSF&, const SSI&) {}
};


} // binary


#define BASE_CLASS(Aux,Class) mpl::if_c<IS_##Aux,binary::Class,binary::EmptyBase<binary::Class> >::type


template<bool IS_EX, bool IS_HA, bool IS_LI>
class BinarySystem
  : public binary::Base,
    public BASE_CLASS(EX,Exact), 
    public BASE_CLASS(HA,Hamiltonian),
    public BASE_CLASS(LI,Liouvillean)
{
public:
  typedef typename BASE_CLASS(EX,Exact)             ExactBase;
  typedef typename BASE_CLASS(HA,Hamiltonian) HamiltonianBase;
  typedef typename BASE_CLASS(LI,Liouvillean) LiouvilleanBase;
  
  typedef structure::Interaction<2> Interaction;

  BinarySystem(Interaction::Ptr);

};


#undef BASE_CLASS


#endif // ELEMENTS_COMPOSITES_BINARYSYSTEM_H_INCLUDED
