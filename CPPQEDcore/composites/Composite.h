// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED

#include "CompositeFwd.h"

#include "Act.h"

#include "BlitzArraySliceIterator.h"
// This is included at this point mainly to pull in necessary TMP tools

#include "SubSystem.h"

#include "details_TMP_helpers.h"

#include <boost/fusion/container/generation/make_list.hpp>

#include <array>

namespace composite {

using boost::fusion::make_list;

using ::size_t;

namespace result_of {

using boost::fusion::result_of::make_list;

} // result_of


template<typename VA>
struct MaxRank : MaxMF<typename MaxMF<VA,SeqLess<mpl::_1,mpl::_2> >::type>::type::type {};



template<int N_RANK>
// Factoring out code that depends only on RANK:
class RankedBase : public structure::QuantumSystem<N_RANK>
{
public:
  static const int RANK=N_RANK;

  typedef boost::shared_ptr<const RankedBase<RANK> > Ptr;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef structure::QuantumSystem<RANK> QS_Base;

  typedef typename QS_Base::Dimensions Dimensions;

  typedef tmptools::Ordinals<RANK> Ordinals;

private:
  // Constructor helper
  static auto fillDimensions(const Frees&);

protected:
  explicit RankedBase(const Frees& frees)
    : QS_Base(fillDimensions(frees)), frees_(frees) {}
  
  const Frees& getFrees() const {return frees_;}

private:
  const Frees frees_;
 
};


template<typename VA>
// VA should model a fusion sequence of Acts
class Base
  : public RankedBase<MaxRank<VA>::value+1>,
    public structure::Averaged<MaxRank<VA>::value+1>
{
public:
  typedef boost::shared_ptr<const Base<VA> > Ptr;
  
  // The calculated RANK
  static const int RANK=MaxRank<VA>::value+1;

private:
  // Compile-time sanity check
  static_assert( composite::CheckMeta<RANK,VA>::type::value == true , "Composite not consistent" );

public:
  // Public types
  typedef VA Acts;

  typedef RankedBase<RANK>          RBase  ;
  typedef structure::Averaged<RANK> Av_Base;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename Av_Base::Averages Averages  ;

  typedef typename RBase::   Frees    Frees;
  typedef typename RBase::Ordinals Ordinals;
  
  template<structure::LiouvilleanAveragedTag>
  static std::ostream& displayKeyLA(std::ostream& , size_t&, const Frees&, const VA& acts);

  template<structure::LiouvilleanAveragedTag>
  class DisplayKey;
  
  template<structure::LiouvilleanAveragedTag>
  static size_t nAvrLA(const Frees& frees, const VA& acts);

  template<structure::LiouvilleanAveragedTag>
  class NAvr;
 
  template<structure::LiouvilleanAveragedTag>
  static const Averages averageLA(double t, const LazyDensityOperator& ldo, const Frees& frees, const VA& acts, size_t numberAvr);

  template<structure::LiouvilleanAveragedTag>
  class Average;

protected:
  // Constructor
  explicit Base(const Frees& frees, const VA& acts)
    : RBase(frees), frees_(RBase::getFrees()), acts_(acts) {}
    
  const VA& getActs() const {return acts_;}

private:
  // Implementing QuantumSystem interface

  double         highestFrequency_v(             ) const;
  std::ostream& displayParameters_v(std::ostream&) const; class DisplayParameters;

  // Implementing Av_Base

  std::ostream& displayKey_v(std::ostream& os, size_t& i)              const {return displayKeyLA<structure::LA_Av>(os,i, frees_,acts_         );}
  size_t              nAvr_v()                                         const {return       nAvrLA<structure::LA_Av>(      frees_,acts_         );}
  const Averages   average_v(double t, const LazyDensityOperator& ldo) const {return    averageLA<structure::LA_Av>(t,ldo,frees_,acts_,nAvr_v());}
  
  void          process_v(Averages&)                           const; class Process;
  std::ostream& display_v(const Averages&, std::ostream&, int) const; class Display;

  const Frees& frees_;
  const VA      acts_;

};


// Constructor helper
template<typename VA>
auto fillFrees(const VA& acts);



template<typename VA>
const typename Base<VA>::Ptr doMake(const VA&);



template<typename VA>
class Exact
  : public structure::Exact<MaxRank<VA>::value+1>
{
private:
  static const int RANK=MaxRank<VA>::value+1;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Exact(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  bool applicableInMaster_v(                     ) const; class IsUnitary;
  void actWithU_v(double, StateVectorLow&, double) const; class ActWithU ;

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename VA>
class Hamiltonian
  : public structure::Hamiltonian<MaxRank<VA>::value+1>
{
private:
  static const int RANK=MaxRank<VA>::value+1;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Hamiltonian(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const; class AddContribution;

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename VA>
class Liouvillean
  : public structure::Liouvillean<MaxRank<VA>::value+1>
{
private:
  static const int RANK=MaxRank<VA>::value+1;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;
  typedef typename quantumdata::Types<RANK>::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef tmptools::Ordinals<RANK> Ordinals;

  typedef typename structure::Liouvillean<RANK>::Rates Rates;

protected:
  Liouvillean(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  std::ostream& displayKey_v(std::ostream& os, size_t& i)              const override {return Base<VA>::template displayKeyLA<structure::LA_Li>(os,i, frees_,acts_         );}
  size_t              nAvr_v()                                         const override {return Base<VA>::template       nAvrLA<structure::LA_Li>(      frees_,acts_         );}
  const Rates      average_v(double t, const LazyDensityOperator& ldo) const override {return Base<VA>::template    averageLA<structure::LA_Li>(t,ldo,frees_,acts_,nAvr_v());}

  void actWithJ_v(double, StateVectorLow&, size_t) const override; class ActWithJ;
  void actWithSuperoperator_v(double, const DensityOperatorLow&, DensityOperatorLow&, size_t) const override; class ActWithSuperoperator;

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename>
class EmptyBase
{
public:
  template<typename VA>
  EmptyBase(const std::array<SubSystemFree,MaxRank<VA>::value+1>&, const VA&) {}
  
};


} // composite



#define BASE_class(Aux,Class) mpl::if_c<IS_##Aux,composite::Class<VA>,composite::EmptyBase<composite::Class<VA> > >::type

/// Class representing a full-fledged composite quantum system defined by a network of \link composite::_ interactions\endlink
/**
 * Assume a system composed of harmonic-oscillator modes and particle motional degrees of freedom layed out in the following way – cf. \ref userguide,
 * the layout means which free-system quantum number corresponds to which index of the multi-array:
 * (Free No. 0) mode (1) motional (2) mode (3) motional (4) mode (5) motional (6) mode (7) motional.
 *
 * Assume further an interaction which couples two modes and two motional degrees of freedom.
 * As a physical example, one might think about two resonator modes with orthogonal axes, and a polarizable particle moving in these two dimensions
 * (the rest of the system can be another particle, and two more modes of the same resonators).
 * Then, the particle may absorb a photon from one mode and suffer recoil in this direction, and emit the photon into the other mode,
 * suffering recoil in that direction as well. Assume that the function
 *
 *     void TwoModesParticle2D_Hamiltonian(const StateVector<4>&, StateVector<4>&);
 *
 * implements this Hamiltonian. In the composite system, this is how we can act with this Hamiltonian between e.g. the indices (4) (0) (7) (1):
 *
 *     void Hamiltonian(const StateVector<8>& psi, StateVector<8>& dpsidt)
 *     {
 *       for_each(BASI_Range<Vector<4,0,7,1> >(psi),
 *                begin<Vector<4,0,7,1> >(dpsidt),
 *                TwoModesParticle2D_Hamiltonian);
 *     }
 *
 * \note the names in this code snippet are approximative
 *
 * It is the task of Composite to systematically keep track of the subsystems and perform the calculations on the necessary slices.
 *
 * Consistency checks for Composite at compile time come on three levels:
 *
 * 1. Composite itself only checks whether all the quantum numbers are addressed by a composite::_ instant. So this is e.g. not allowed:
 *
 *        composite::result_of::Make<composite::_<0,1>,composite::_<0,2>,composite::_<0,3>,composite::_<3,2,5,1> >::type
 *
 *    because 4 is not addressed and hence the Composite object has no way to figure out what kind of structure::Free object is there.
 *
 * 2. The instantiated \link blitzplusplus::basi::Iterator slice iterators\endlink do some further checks for each composite::_ instant individually:
 *   - the total arity must be greater than or equal to the arity of composite::_
 *   - composite::_ must not “contain” duplicated “elements” (e.g. composite::_<3,2,3,1> not allowed)
 *   - each element in composite::_ must be smaller than the total arity
 *
 * 3. Finally, tmptools::Vector checks for the non-negativity of each element (as it is supposed to be a non-negative compile-time vector).
 *
 * \note The last condition incidentally makes that the three checks performed under 2. by the \link blitzplusplus::basi::Iterator slice interators\endlink
 * for each composite::_ are redundant (the first condition follows from the last two plus the nonnegativity of tmptools::Vector).
 * However, these checks are still necessary for the slice iterators
 * because they accept classes other than tmptools::Vector for the specification of retained index positions.
 *
 * This is followed by a check at runtime, when the actual elements become available, whether the legs of the composite::_ objects are consistent among each other.
 * Cf. composite::FillFrees::Inner::operator().
 *
 * \tparam VA should model a \refBoost{Boost.Fusion list,fusion/doc/html/fusion/container/list.html} of composite::_ objects
 * \tparam IS_EX governs whether the class should inherit from composite::Exact
 * \tparam IS_HA governs whether the class should inherit from composite::Hamiltonian
 * \tparam IS_LI governs whether the class should inherit from composite::Liouvillean
 */
template<typename VA, bool IS_EX=true, bool IS_HA=true, bool IS_LI=true>
// VA should model a fusion sequence of Acts
class Composite
  : public composite::Base<VA>,
    public BASE_class(EX,Exact),
    public BASE_class(HA,Hamiltonian),
    public BASE_class(LI,Liouvillean)
{
public:
  typedef composite::Base<VA> Base;
  
  typedef typename BASE_class(EX,Exact)             ExactBase;
  typedef typename BASE_class(HA,Hamiltonian) HamiltonianBase;
  typedef typename BASE_class(LI,Liouvillean) LiouvilleanBase;
  
  typedef typename composite::Base<VA>::Frees Frees;
  
  // The calculated RANK
  static const int RANK=Base::RANK;

private:
  using Base::getFrees; using Base::getActs ;
  
public:
  Composite(const Frees& frees, const VA& acts)
    : Base(frees ,acts),
      ExactBase      (getFrees(),getActs()),
      HamiltonianBase(getFrees(),getActs()),
      LiouvilleanBase(getFrees(),getActs()) {}

  // Constructor
  explicit Composite(const VA& acts)
    : Composite(composite::fillFrees(acts),acts) {}

  friend const typename composite::Base<VA>::Ptr composite::doMake<VA>(const VA&);

  // Note that Frees and Acts are stored by value in Base

};

#undef BASE_class


// The following provides a much more convenient interface:

namespace composite {

namespace result_of {


template<bool IS_EX, bool IS_HA, bool IS_LI, typename... Acts>
struct MakeConcrete : boost::mpl::identity<Composite<typename make_list<Acts...>::type, IS_EX, IS_HA, IS_LI> > {};


template<typename... Acts>
struct Make : boost::mpl::identity<typename Base<typename make_list<Acts...>::type>::Ptr> {};



} // result_of


template<typename... Acts>
auto make(const Acts&... acts)
{
  return doMake(make_list(acts...));
}


} // composite


#endif // CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED
