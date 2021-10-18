// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_COMPOSITES_BINARYSYSTEM_H_INCLUDED
#define CPPQEDCORE_COMPOSITES_BINARYSYSTEM_H_INCLUDED

#include "LazyDensityOperator.h"
#include "QuantumSystem.h"
#include "SubSystem.h"

namespace mpl=boost::mpl;

/// Auxiliary tools for BinarySystem
namespace binary {

typedef ::structure::Interaction<2> Interaction; ///< Binary interaction
using InteractionPtr = ::structure::InteractionPtr<2>;

using StateVectorLow = ::structure::StateVectorLow<2>;
using DensityOperatorLow = ::structure::DensityOperatorLow<2>;

using LazyDensityOperator = ::quantumdata::LazyDensityOperator<2>;

typedef composite::SubSystemFree SSF; ///< Convenience typedef
typedef composite::SubSystemsInteraction<2> SSI; ///< Convenience typedef

using structure::Averages; using structure::Rates;

/// Outfactored common functionality of Liouvillean and Averaged
template<structure::LiouvilleanAveragedTag>
std::ostream& streamKey(std::ostream&, size_t&, const SSF& free0, const SSF& free1, const SSI& ia);

/// Outfactored common functionality of Liouvillean and Averaged
template<structure::LiouvilleanAveragedTag>
size_t nAvr(const SSF& free0, const SSF& free1, const SSI& ia);

/// Outfactored common functionality of Liouvillean and Averaged
template<structure::LiouvilleanAveragedTag>
const Averages
average(double t, const quantumdata::LazyDensityOperator<2>& ldo, const SSF& free0, const SSF& free1, const SSI& ia, size_t numberAvr);


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
  typedef structure::Averaged<1> Av1;
  typedef structure::Averaged<2> Av2;

  /// Constructor from an #Interaction instant
  explicit Base(InteractionPtr);

  /// \name Getters
  //@{
  const SSF& getFree0() const {return free0_;}
  const SSF& getFree1() const {return free1_;}
  const SSI& getIA () const {return ia_;}
  //@}

private:
  double highestFrequency_v() const;

  std::ostream& streamParameters_v(std::ostream&) const;

  size_t nAvr_v() const {return binary::nAvr <structure::LA_Av>(free0_,free1_,ia_);}
  
  const Averages average_v(double t, const LazyDensityOperator& ldo) const {return binary::average <structure::LA_Av>(t,ldo,free0_,free1_,ia_,nAvr());}
  
  void process_v(Averages&) const;
  
  std::ostream& stream_v(const Averages&, std::ostream&, int) const;
  
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const {return binary::streamKey<structure::LA_Av>(os,i, free0_,free1_,ia_);}

  const SSF free0_, free1_;

  const SSI ia_;
  
};

typedef std::shared_ptr<const Base> Ptr; ///< Convenience typedef

/// Maker function for BinarySystem
/**
 * Uses runtime dispatching to select the suitable class-composed BinarySystem on the basis of the characteristics of the components
 * (the two \link structure::Free free systems\endlink and the #Interaction): whether they derive from structure::Exact, structure::Hamiltonian, structure::Liouvillean.
 * If any of the components derive from structure::Exact, then the whole BinarySystem has to derive from structure::Exact, and so on.
 */
const Ptr make(InteractionPtr);


#define CLASS_HEADER(Class) class Class : public structure::Class<2>

#define CLASS_BODY_PART(Class,Aux) public:                                   \
  typedef structure::Class<1> Aux##1;                                        \
  typedef structure::Class<2> Aux##2;                                        \
                                                                             \
  Class(const SSF& free0, const SSF& free1, const SSI& ia) : free0_(free0), free1_(free1), ia_(ia) {} \
                                                                             \
private:                                                                     \
  const SSF &free0_, &free1_;                                                \
  const SSI &ia_;                                                            \


/// Implements the structure::Exact interface for a BinarySystem along the same lines as Base implements the structure::Averaged interface
CLASS_HEADER(Exact)
{
  CLASS_BODY_PART(Exact,Ex)

  bool applicableInMaster_v() const;

  void actWithU_v(double, StateVectorLow&, double) const;

};


/// Implements the structure::Hamiltonian interface for a BinarySystem
CLASS_HEADER(Hamiltonian)
{
  CLASS_BODY_PART(Hamiltonian,Ha)

  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const;

};


/// Implements the structure::Liouvillean interface for a BinarySystem
CLASS_HEADER(Liouvillean)
{
  CLASS_BODY_PART(Liouvillean,Li)

  void actWithJ_v(double, StateVectorLow&, size_t) const override;
  
  void actWithSuperoperator_v(double, const DensityOperatorLow&, DensityOperatorLow&, size_t) const override;

  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return binary::streamKey<structure::LA_Li>(os,i, free0_,free1_,ia_);}
  
  size_t nAvr_v() const override {return binary::nAvr<structure::LA_Li>(free0_,free1_,ia_);}
  
  const Rates average_v(double t, const LazyDensityOperator& ldo) const override {return binary::average<structure::LA_Li>(t,ldo,free0_,free1_,ia_,nAvr());}

};


#undef CLASS_BODY_PART
#undef CLASS_HEADER

/// Helper for class composition of BinarySystem
template<typename>
class EmptyBase
{
public:
  EmptyBase(const SSF&, const SSF&, const SSI&) {}
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
