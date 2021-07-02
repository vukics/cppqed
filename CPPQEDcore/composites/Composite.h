// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED

#include "SubSystem.h"

#include "SliceIterator.h"

#include <array>


namespace composite {


using ::size_t; using ::structure::Averages; using ::structure::Rates;


/// Helper class to Composite (an Act)
/**
 * It combines a set of \link retainedindexpositionsdefined retained index positions\endlink (via its tmptools::Vector base) with a certain structure::Interaction element
 * (via the SubSystemsInteraction base) of the corresponding arity (equalling the number of retained index positions) – this means that the given structure::Interaction
 * element will act on the given retained index positions.
 *
 * Composite expects a \refBoost{Boost.Fusion,fusion} sequence as its template argument `Acts`.
 *
 * \tparam RetainedAxes has the same role as in tmptools::Vector
 */
template<int... RetainedAxes>
class _
  : public SubSystemsInteraction<sizeof...(RetainedAxes)>
{
public:
  explicit _(typename structure::InteractionPtr<sizeof...(RetainedAxes)> ia) : SubSystemsInteraction<sizeof...(RetainedAxes)>(ia) {}

  constexpr operator tmptools::Vector<RetainedAxes...>() const {return tmptools::vector<RetainedAxes...>;};
  
};


template<typename Acts>
constexpr int calculateRank(Acts acts = Acts{})
{
  return hana::maximum(hana::maximum(hana::maximum(acts, [&](auto&& a, auto&& b) {return hana::maximum(a)<hana::maximum(b);}))) + 1;
}


template<int RANK>
// Factoring out code that depends only on RANK:
class RankedBase : public structure::QuantumSystem<RANK>
{
public:
  static const int N_RANK = RANK;

  using Frees = std::array<SubSystemFree,RANK> ;

protected:
  explicit RankedBase(const Frees& frees)
    : structure::QuantumSystem<RANK>([&] () {
        typename structure::QuantumSystem<RANK>::Dimensions res;
        hana::for_each(tmptools::ordinals<RANK>,[&](auto t) {res(t)=frees[t].get()->getDimension();});
        return res;
    }()), frees_(frees) {}
  
  const Frees& getFrees() const {return frees_;}

private:
  const Frees frees_;
 
};


template<typename Acts> // Acts should be a hana::tuple of Acts
class Base
  : public RankedBase<calculateRank<Acts>()>,
    public structure::Averaged<calculateRank<Acts>()>
{
public:
  // The calculated RANK
  static constexpr int RANK = calculateRank<Acts>;

private:
  using Frees = typename RankedBase<RANK>::Frees; using RankedBase<RANK>::getFrees;
  
  // Compile-time sanity check: each ordinal up to RANK has to be contained by at least one act.
  // (That is, no free system can be left out of the network of interactions)
  static_assert( hana::fold (tmptools::ordinals<RANK>, hana::true_c, [](auto state, auto ordinal) {
    return state && hana::fold (Acts{}, hana::false_c, [o=ordinal](auto state, auto act) {return state || hana::contains(act,o);}); }
  ) == true , "Composite not consistent" );

public:
  template<structure::LiouvilleanAveragedTag>
  static std::ostream& streamKeyLA(std::ostream&, size_t&, const Frees&, const Acts& acts);

  template<structure::LiouvilleanAveragedTag>
  static size_t nAvrLA(const Frees& frees, const Acts& acts);

  template<structure::LiouvilleanAveragedTag>
  static const structure::Averages averageLA(double t, const quantumdata::LazyDensityOperator<RANK>& ldo, const Frees& frees, const Acts& acts, size_t numberAvr);

protected:
  // Constructor
  explicit Base(const Frees& frees, const Acts& acts) : RankedBase<RANK>(frees), acts_(acts) {}
    
  const Acts& getActs() const {return acts_;}

private:
  // Implementing QuantumSystem interface

  double highestFrequency_v( ) const;
  std::ostream& streamParameters_v(std::ostream&) const;

  // Implementing Av_Base

  std::ostream& streamKey_v(std::ostream& os, size_t& i) const {return streamKeyLA<structure::LA_Av>(os,i, getFrees() ,acts_ );}
  size_t nAvr_v() const {return nAvrLA<structure::LA_Av>( getFrees(),acts_ );}
  const structure::Averages average_v(double t, const quantumdata::LazyDensityOperator<RANK>& ldo) const override {return averageLA<structure::LA_Av>(t,ldo,getFrees(),acts_,nAvr_v());}
  
  void process_v(structure::Averages&) const;
  std::ostream& stream_v(const structure::Averages&, std::ostream&, int) const;

  const Acts acts_;

};


// Constructor helper
template<typename Acts>
auto fillFrees(const Acts& acts);



template<typename Acts>
const typename Base<Acts>::Ptr doMake(const Acts&);



template<typename Acts>
class Exact
  : public structure::Exact<MaxRank_v<Acts>+1>
{
private:
  static const int RANK=MaxRank_v<Acts>+1;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef quantumdata::StateVectorLow<RANK> StateVectorLow;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Exact(const Frees& frees, const Acts& acts) : frees_(frees), acts_(acts) {}

private:
  void actWithU_v(double, StateVectorLow&, double) const;
  bool applicableInMaster_v( ) const;

  const Frees& frees_;
  const Acts & acts_;

};


template<typename Acts>
class Hamiltonian
  : public structure::Hamiltonian<MaxRank_v<Acts>+1>
{
private:
  static const int RANK=MaxRank_v<Acts>+1;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef quantumdata::StateVectorLow<RANK> StateVectorLow;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Hamiltonian(const Frees& frees, const Acts& acts) : frees_(frees), acts_(acts) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const;

  const Frees& frees_;
  const Acts & acts_;

};


template<typename Acts>
class Liouvillean
  : public structure::Liouvillean<MaxRank_v<Acts>+1>
{
private:
  static const int RANK=MaxRank_v<Acts>+1;

  typedef std::array<SubSystemFree,RANK> Frees;

  typedef quantumdata::StateVectorLow<RANK> StateVectorLow;
  typedef quantumdata::DensityOperatorLow<RANK> DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Liouvillean(const Frees& frees, const Acts& acts) : frees_(frees), acts_(acts) {}

private:
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return Base<Acts>::template streamKeyLA<structure::LA_Li>(os,i, frees_,acts_);}
  
  size_t nAvr_v() const override {return Base<Acts>::template nAvrLA<structure::LA_Li>(frees_,acts_);}
  
  const Rates average_v(double t, const LazyDensityOperator& ldo) const override {return Base<Acts>::template averageLA<structure::LA_Li>(t,ldo,frees_,acts_,nAvr_v());}

  void actWithJ_v(double, StateVectorLow&, size_t) const override;
  void actWithSuperoperator_v(double, const DensityOperatorLow&, DensityOperatorLow&, size_t) const override;

  const Frees& frees_;
  const Acts & acts_;

};


template<typename>
class EmptyBase
{
public:
  template<typename Acts>
  EmptyBase(const std::array<SubSystemFree,MaxRank_v<Acts>+1>&, const Acts&) {}
  
};


} // composite



#define BASE_class(Aux,Class) std::conditional_t<IS_##Aux,composite::Class<Acts>,composite::EmptyBase<composite::Class<Acts>>>

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
 * \tparam Acts should model a \refBoost{Boost.Fusion list,fusion/doc/html/fusion/container/list.html} of composite::_ objects
 * \tparam IS_EX governs whether the class should inherit from composite::Exact
 * \tparam IS_HA governs whether the class should inherit from composite::Hamiltonian
 * \tparam IS_LI governs whether the class should inherit from composite::Liouvillean
 */
template<typename Acts, bool IS_EX=true, bool IS_HA=true, bool IS_LI=true>
// Acts should model a fusion sequence of Acts
class Composite
  : public composite::Base<Acts>,
    public BASE_class(EX,Exact),
    public BASE_class(HA,Hamiltonian),
    public BASE_class(LI,Liouvillean)
{
public:
  typedef composite::Base<Acts> Base;
  
  typedef BASE_class(EX,Exact) ExactBase;
  typedef BASE_class(HA,Hamiltonian) HamiltonianBase;
  typedef BASE_class(LI,Liouvillean) LiouvilleanBase;
  
  typedef typename composite::Base<Acts>::Frees Frees;
  
  // The calculated RANK
  static const int RANK=Base::RANK;

private:
  using Base::getFrees; using Base::getActs ;
  
public:
  Composite(const Frees& frees, const Acts& acts)
    : Base(frees ,acts),
      ExactBase(getFrees(),getActs()),
      HamiltonianBase(getFrees(),getActs()),
      LiouvilleanBase(getFrees(),getActs()) {}

  // Constructor
  explicit Composite(const Acts& acts)
    : Composite(composite::fillFrees(acts),acts) {}

  friend const typename composite::Base<Acts>::Ptr composite::doMake<Acts>(const Acts&);

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


/// Template alias for backward compatibility
template<int... RetainedAxes>
using Act = composite::_<RetainedAxes...>;


#endif // CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED
