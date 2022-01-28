// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_COMPOSITES_BINARYSYSTEM_H_INCLUDED
#define CPPQEDCORE_COMPOSITES_BINARYSYSTEM_H_INCLUDED

#include "LazyDensityOperator.h"
#include "Structure.h"

#include "Interaction.h"

namespace mpl=boost::mpl;

/// Auxiliary tools for BinarySystem
namespace binary {

typedef tmptools::Vector<0> V0;
typedef tmptools::Vector<1> V1;

using structure::FreePtr;

typedef ::structure::Interaction<2> Interaction; ///< Binary interaction
using InteractionPtr = ::structure::InteractionPtr<2>;

using StateVectorLow = ::structure::StateVectorLow<2>;
using DensityOperatorLow = ::structure::DensityOperatorLow<2>;

using LazyDensityOperator = ::quantumdata::LazyDensityOperator<2>;


using structure::Averages; using structure::Rates;

/// Outfactored common functionality of Liouvillean and Averaged
//@{
template<structure::LiouvilleanAveragedTag LA>
std::ostream& streamKey(InteractionPtr ia, std::ostream& os, size_t& i)
{  
  return ::structure::streamKey<LA>(ia,
                                    ::structure::streamKey<LA>((*ia)[1],
                                                               ::structure::streamKey<LA>((*ia)[0],
                                                                                          os<<"Binary system\n",i),i),i);
}


template<structure::LiouvilleanAveragedTag LA>
size_t nAvr(InteractionPtr ia)
{
  return structure::nAvr<LA>((*ia)[0]) + structure::nAvr<LA>((*ia)[1]) + structure::nAvr<LA>(ia);
}


template<structure::LiouvilleanAveragedTag LA>
const Averages average(InteractionPtr ia, double t, const quantumdata::LazyDensityOperator<2>& ldo, size_t numberAvr)
{
  using boost::copy;

  const Averages
    a0 {quantumdata::partialTrace<V0>(ldo,[&](const auto& m){return ::structure::average<LA>((*ia)[0],t,m);})},
    a1 {quantumdata::partialTrace<V1>(ldo,[&](const auto& m){return ::structure::average<LA>((*ia)[1],t,m);})},
    a01{::structure::average<LA>(ia,t,ldo)};

  Averages a(numberAvr);

  copy(a01,copy(a1,copy(a0,a.begin())));

  return a;
}

//@}

/// Common base for all class-composed BinarySystem%s
/**
 * It implements only the structure::Averaged and structure::QuantumSystem interfaces, since all the other characteristics
 * (structure::Exact, structure::Hamiltonian, structure::Liouvillean) will be added by class composition.
 *
 * \see explanation @ make()
 *
 * \par Example
 *
 * structure::Averaged::average is implemented by taking the corresponding LazyDensityOperator unary slices for the `average` functions of the two frees
 * (if they derive from structure::Average at all), and passing the complete binary LazyDensityOperator to the `average` function of the interaction component.
 *
 * \see quantumdata::partialTrace
 *
 */
class Base
  : public structure::QuantumSystem<2>,
    public structure::Averaged <2>
{
protected:
  explicit Base(InteractionPtr ia) : QuantumSystem<2>(ia->getDimensions()), ia_(ia) {}

private:
  double highestFrequency_v() const override;

  std::ostream& streamParameters_v(std::ostream&) const override;

  size_t nAvr_v() const override {return binary::nAvr <structure::LA_Av>(ia_);}
  
  const Averages average_v(double t, const LazyDensityOperator& ldo) const override {return binary::average<structure::LA_Av>(ia_,t,ldo,nAvr());}
  
  void process_v(Averages&) const override;
  
  std::ostream& stream_v(const Averages&, std::ostream&, int) const override;
  
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return binary::streamKey<structure::LA_Av>(ia_,os,i);}

protected:
  const InteractionPtr ia_;
  
};

typedef std::shared_ptr<const Base> Ptr; ///< Convenience typedef

/// Maker function for BinarySystem
/**
 * Uses runtime dispatching to select the suitable class-composed BinarySystem on the basis of the characteristics of the components
 * (the two \link structure::Free free systems\endlink and the #Interaction): whether they derive from structure::Exact, structure::Hamiltonian, structure::Liouvillean.
 * If any of the components derive from structure::Exact, then the whole BinarySystem has to derive from structure::Exact, and so on.
 */
const Ptr make(InteractionPtr);


/// Implements the structure::Exact interface for a BinarySystem along the same lines as Base implements the structure::Averaged interface
class Exact : public structure::Exact<2>
{
public:
  explicit Exact(InteractionPtr ia) : ia_{ia} {}
  
private:
  bool applicableInMaster_v() const override;

  void actWithU_v(double, StateVectorLow&, double) const override;
  
  const InteractionPtr ia_;

};


/// Implements the structure::Hamiltonian interface for a BinarySystem
class Hamiltonian : public structure::Hamiltonian<2>
{
public:
  explicit Hamiltonian(InteractionPtr ia) : ia_{ia} {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const override;

  const InteractionPtr ia_;
  
};


/// Implements the structure::Liouvillean interface for a BinarySystem
class Liouvillean : public structure::Liouvillean<2>
{
public:
  explicit Liouvillean(InteractionPtr ia) : ia_{ia} {}

private:
  void actWithJ_v(double, StateVectorLow&, size_t) const override;
  
  void actWithSuperoperator_v(double, const DensityOperatorLow&, DensityOperatorLow&, size_t) const override;

  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return binary::streamKey<structure::LA_Li>(ia_,os,i);}
  
  size_t nAvr_v() const override {return binary::nAvr<structure::LA_Li>(ia_);}
  
  const Rates average_v(double t, const LazyDensityOperator& ldo) const override {return binary::average<structure::LA_Li>(ia_,t,ldo,nAvr());}

  const InteractionPtr ia_;

};


/// Helper for class composition of BinarySystem
template<typename>
class EmptyBase
{
public:
  EmptyBase(InteractionPtr) {}
};


} // binary


#define BASE_class(Aux,Class) std::conditional_t<IS_##Aux,binary::Class,binary::EmptyBase<binary::Class> >


/// Implements the simplest composite system: a binary where a single binary::Interaction couples two \link structure::Free free systems\endlink
/**
 * The class is meant as a substitute for the full Composite in this simplest case, for saving compile-time (and, maybe, also runtime) resources.
 *
 * It inherits unconditionally from binary::Base, so that it has the structure::QuantumSystem interface necessary to be evolved by
 * \link #quantumtrajectory quantum trajectories\endlink. All other \link #structure interfaces\endlink are added conditionally:
 *
 * \tparam IS_EX governs whether the class should inherit from binary::Exact
 * \tparam IS_HA governs whether the class should inherit from binary::Hamiltonian
 * \tparam IS_LI governs whether the class should inherit from binary::Liouvillean
 *
 * If `IS_EX` is `true`, the class inherits from binary::Exact, otherwise from binary::EmptyBase, and so on.
 *
 */
template<bool IS_EX, bool IS_HA, bool IS_LI>
class BinarySystem
  : public binary::Base,
    public BASE_class(EX,Exact), 
    public BASE_class(HA,Hamiltonian),
    public BASE_class(LI,Liouvillean)
{
public:
  typedef BASE_class(EX,Exact) ExactBase;
  
  typedef BASE_class(HA,Hamiltonian) HamiltonianBase;
  
  typedef BASE_class(LI,Liouvillean) LiouvilleanBase;
  
  explicit BinarySystem(binary::InteractionPtr);

};


#undef BASE_class


#endif // CPPQEDCORE_COMPOSITES_BINARYSYSTEM_H_INCLUDED
