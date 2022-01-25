// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED

#include "SubSystem.h"

#include "Algorithm.h"
#include "SliceIterator.h"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include <array>



namespace composite {


struct ConsistencyException : std::logic_error
{
  ConsistencyException(int idxx, int ii) : std::logic_error(std::to_string(idxx)+" "+std::to_string(ii)), idx(idxx), i(ii) {}

  const int idx, i;
};


using ::size_t; using ::structure::Averages; using ::structure::Rates;

template<int RANK>
using Frees = std::array<SubSystemFree,RANK> ;


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
  explicit _(structure::InteractionPtr<sizeof...(RetainedAxes)> ia) : SubSystemsInteraction<sizeof...(RetainedAxes)>(ia) {}

  static constexpr auto retainedAxes=tmptools::vector<RetainedAxes...>;
  
};


template<typename... Acts>
constexpr auto rank_v=hana::maximum(hana::flatten(hana::make_tuple(Acts::retainedAxes...)))+1;


template<int RANK>
// Factoring out code that depends only on RANK:
class RankedBase : public structure::QuantumSystem<RANK>
{
public:
  static const int N_RANK = RANK;

protected:
  explicit RankedBase(const Frees<RANK>& frees)
    : structure::QuantumSystem<RANK>([&] () { // calculate the dimensions
        typename structure::QuantumSystem<RANK>::Dimensions res;
        hana::for_each(tmptools::ordinals<RANK>,[&](auto t) {res(t)=frees[t]->getDimension();});
        return res;
    }()), frees_(frees) {}
  
  const auto& getFrees() const {return frees_;}

private:
  const Frees<RANK> frees_;
 
};


template<typename... Acts> // Acts should be a hana::tuple of Acts
class Base
  : public RankedBase<rank_v<Acts...>>,
    public structure::Averaged<rank_v<Acts...>>
{
public:
  // The calculated RANK
  static constexpr int RANK = rank_v<Acts...>;

protected:
  using RankedBase<RANK>::getFrees;
  
  // Compile-time sanity check: each ordinal up to RANK has to be contained by at least one act.
  // (That is, no free system can be left out of the network of interactions)
  static_assert( hana::fold (tmptools::ordinals<RANK>, hana::true_c, [](auto state, auto ordinal) {
    return state && hana::contains(hana::flatten(hana::make_tuple(Acts::retainedAxes...)),ordinal);
  }) == true , "Composite not consistent" );

public:
  template<structure::LiouvilleanAveragedTag LA>
  static std::ostream& streamKeyLA(std::ostream& os, size_t& i, const Frees<RANK>& frees, const hana::tuple<Acts...>& acts)
  {
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {frees[idx].template streamKey<LA>(os,i);});
    hana::for_each(acts,[&](const auto& act) {act.template streamKey<LA>(os,i);});    
    return os;
  }
  
  
  template<structure::LiouvilleanAveragedTag LA>
  static size_t nAvrLA(const Frees<RANK>& frees, const hana::tuple<Acts...>& acts)
  {
    
    return hana::fold(acts,
                      hana::fold(tmptools::ordinals<RANK>,0,[&](size_t res, auto idx) {return res + frees[idx].template nAvr<LA>();}),
                      [&](size_t s, const auto& act) {return s+act.template nAvr<LA>();});
  }
  

#ifndef NDEBUG
#pragma GCC warning "TODO: This solution is a bit insane, could be solved more effectively with blitz::Array slice"
#endif // NDEBUG
  template<structure::LiouvilleanAveragedTag LA>
  static const structure::Averages averageLA(double t, const quantumdata::LazyDensityOperator<RANK>& ldo,
                                             const Frees<RANK>& frees, const hana::tuple<Acts...>& acts, size_t numberAvr)
  {
    std::list<Averages> seqAverages{RANK+hana::size(acts)}; // Averages res(numberAvr); size_t resIdx=0;
    
    {
      typename std::list<Averages>::iterator iter(seqAverages.begin());
      
      const auto lambda=[&](auto av, auto v) {
        iter++->reference(quantumdata::partialTrace<std::decay_t<decltype(v)>>(ldo,[&](const auto& ldoS) {
          return structure::average(av,t,ldoS);
        }));
      };
      
      auto tag=structure::LiouvilleanAveragedTag_<LA>();
      
      hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees[idx].getLA(tag),tmptools::vector<idx>);});
      
      hana::for_each(acts,[&](const auto& act) {lambda(act.getLA(tag),act.retainedAxes);});
      
    }
    
    Averages res(numberAvr); res=0;
    return cppqedutils::concatenate(seqAverages,res);
    
  }


protected:
  // Constructor
  explicit Base(const Frees<RANK>& frees, const Acts&... acts) : RankedBase<RANK>{frees}, acts_{acts...} {}
    
  const auto& getActs() const {return acts_;}

private:
  // Implementing QuantumSystem interface

  double highestFrequency_v( ) const override
  {
    return std::max(
      getFrees()[hana::maximum(tmptools::ordinals<RANK>,[this](auto idx1, auto idx2) {
        return  getFrees()[idx1]->highestFrequency() < getFrees()[idx2]->highestFrequency();
      })]->highestFrequency()
      ,
      hana::maximum(acts_,[&](const auto& act1, const auto& act2) {return act1->highestFrequency() < act2->highestFrequency();})->highestFrequency()
    );
  }
  
  std::ostream& streamParameters_v(std::ostream& os) const override
  {
    os<<"Composite\nDimensions: "<<RankedBase<RANK>::getDimensions()<<". Total: "<<RankedBase<RANK>::getTotalDimension()<<std::endl;
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {
      os<<"Subsystem Nr. "<<idx<<std::endl;
      getFrees()[idx]->streamParameters(os);
    });
    
    hana::for_each(acts_, [&](const auto& act) {
      hana::for_each(act.retainedAxes,[&](auto idx) {os<<idx<<" - ";});
      os<<"Interaction\n";
      act->streamParameters(os);
    });
    
    return os;
  }
  

  // Implementing the Averaged base

  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return streamKeyLA<structure::LA_Av>(os,i, getFrees() ,acts_ );}
  size_t nAvr_v() const override {return nAvrLA<structure::LA_Av>( getFrees(),acts_ );}
  const structure::Averages average_v(double t, const quantumdata::LazyDensityOperator<RANK>& ldo) const override {return averageLA<structure::LA_Av>(t,ldo,getFrees(),acts_,nAvr_v());}
  
  void process_v(structure::Averages& avr) const override
  {
    ptrdiff_t l=-1, u=0;

    const auto lambda=[&](auto av) {
      if (av && (u=l+av->nAvr())>l) {
        Averages temp(avr(blitz::Range(l+1,u)));
        av->process(temp);
        l=u;
      }
    };
  
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(getFrees()[idx].getAv());});
    
    hana::for_each(acts_,[&](const auto& act) {lambda(act.getAv());});
  }
  
  
  std::ostream& stream_v(const structure::Averages& avr, std::ostream& os, int precision) const override
  {
    ptrdiff_t l=-1, u=0;
  
    const auto lambda=[&](auto av){
      if (av && (u=l+av->nAvr())>l) {
        av->stream(avr(blitz::Range(l+1,u)),os,precision);
        l=u;
      }
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(getFrees()[idx].getAv());});
    
    hana::for_each(acts_,[&](const auto& act) {lambda(act.getAv());});
    
    return os;
    
  }
  
  
  const hana::tuple<Acts...> acts_;

};




template<typename... Acts>
class Exact
  : public structure::Exact<rank_v<Acts...>>
{
private:
  static const int RANK=rank_v<Acts...>;

protected:
  Exact(const Frees<RANK>& frees, const hana::tuple<Acts...>& acts) : frees_(frees), acts_(acts) {}

private:
  void actWithU_v(double t, quantumdata::StateVectorLow<RANK>& psi, double t0) const override
  {
    const auto lambda=[&](auto ex, auto v) {if (ex) for (auto& psiS : cppqedutils::sliceiterator::fullRange<decltype(v)>(psi)) ex->actWithU(t,psiS,t0);};
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx].getEx(),tmptools::vector<idx>);});
    
    hana::for_each(acts_,[&](const auto& act) {lambda(act.getEx(),act.retainedAxes);});
    
  }
  
  
  bool applicableInMaster_v( ) const override
  {
    if (hana::fold(tmptools::ordinals<RANK>,true,[&](bool s, auto idx) {return s && frees_[idx].applicableInMaster();}))
      return hana::fold(acts_,true,[&](bool s, const auto& act) {return s && act.applicableInMaster();});
    else return false;
  }
  

  const Frees<RANK>& frees_;
  const hana::tuple<Acts...> & acts_;

};


template<typename... Acts>
class Hamiltonian
  : public structure::Hamiltonian<rank_v<Acts...>>
{
private:
  static const int RANK=rank_v<Acts...>;

protected:
  Hamiltonian(const Frees<RANK>& frees, const hana::tuple<Acts...>& acts) : frees_(frees), acts_(acts) {}

private:
  void addContribution_v(double t, const quantumdata::StateVectorLow<RANK>& psi, quantumdata::StateVectorLow<RANK>& dpsidt, double t0) const override
  {
    const auto lambda=[&](auto ha, auto v) {
      if (ha)
        boost::for_each(cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(psi),
                        cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(dpsidt),
                        [&](const auto& psiS, auto& dpsidtS) {ha->addContribution(t,psiS,dpsidtS,t0);});
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx].getHa(),tmptools::vector<idx>);});

    hana::for_each(acts_,[&](const auto& act) {lambda(act.getHa(),act.retainedAxes);});
    
  }
  
  const Frees<RANK>& frees_;
  const hana::tuple<Acts...> & acts_;

};


template<typename... Acts>
class Liouvillean
  : public structure::Liouvillean<rank_v<Acts...>>
{
private:
  static const int RANK=rank_v<Acts...>;

protected:
  Liouvillean(const Frees<RANK>& frees, const hana::tuple<Acts...>& acts) : frees_(frees), acts_(acts) {}

private:
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return Base<Acts...>::template streamKeyLA<structure::LA_Li>(os,i, frees_,acts_);}
  
  size_t nAvr_v() const override {return Base<Acts...>::template nAvrLA<structure::LA_Li>(frees_,acts_);}
  
  const Rates average_v(double t, const quantumdata::LazyDensityOperator<RANK>& ldo) const override {return Base<Acts...>::template averageLA<structure::LA_Li>(t,ldo,frees_,acts_,nAvr_v());}

  void actWithJ_v(double t, quantumdata::StateVectorLow<RANK>& psi, size_t ordoJump) const override
  {
    bool flag=false;
    
    const auto lambda=[&](auto li, auto v) {
      if (!flag && li) {
        size_t n=li->nAvr();
        if (ordoJump<n) {
          for (auto& psiS : cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(psi)) li->actWithJ(t,psiS,ordoJump);
          flag=true;
        }
        ordoJump-=n;  
      }
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx].getLi(),tmptools::vector<idx>);});
    
    hana::for_each(acts_,[&](const auto& act) {lambda(act.getLi(),act.retainedAxes);});
    
  }
  
  
  void actWithSuperoperator_v(double t, const quantumdata::DensityOperatorLow<RANK>& rho, quantumdata::DensityOperatorLow<RANK>& drhodt, size_t ordoJump) const override
  {
    bool flag=false;
    
    const auto lambda=[&](auto li, auto v) {
      if (!flag && li) {
        size_t n=li->nAvr();
        if (ordoJump<n) {
          using ExtendedV=decltype(hana::concat(v, /* hana::transform(v,[&](auto element) constexpr {return element+RANK;})*/ v+RANK ));

          boost::for_each(cppqedutils::sliceiterator::fullRange<ExtendedV>(rho),cppqedutils::sliceiterator::fullRange<ExtendedV>(drhodt),
                          [&](const auto& rhoS, auto& drhodtS) {li->actWithSuperoperator(t,rhoS,drhodtS,ordoJump);});
          flag=true;
        }
        ordoJump-=n;  
      }
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx].getLi(),tmptools::vector<idx>);});
    
    hana::for_each(acts_,[&](const auto& act) {lambda(act.getLi(),act.retainedAxes);});
    
  }
  
  const Frees<RANK>& frees_;
  const hana::tuple<Acts...> & acts_;

};


template<typename>
class EmptyBase
{
public:
  template<typename Dummy, typename Acts>
  EmptyBase(const Dummy&, const Acts&) {}
  
};


} // composite



#define BASE_class(Aux,Class) std::conditional_t<IS_##Aux,composite::Class<Acts...>,composite::EmptyBase<composite::Class<Acts...>>>


template<bool IS_EX, bool IS_HA, bool IS_LI, typename... Acts>
// Acts should model a fusion sequence of Acts
class Composite
  : public composite::Base<Acts...>,
    public BASE_class(EX,Exact),
    public BASE_class(HA,Hamiltonian),
    public BASE_class(LI,Liouvillean)
{
public:
  typedef composite::Base<Acts...> Base;
  
  typedef BASE_class(EX,Exact) ExactBase;
  typedef BASE_class(HA,Hamiltonian) HamiltonianBase;
  typedef BASE_class(LI,Liouvillean) LiouvilleanBase;
  
  // The calculated RANK
  static const int RANK=Base::RANK;

private:
  using Base::getFrees; using Base::getActs ;
  
public:
  // Constructor
  Composite(const composite::Frees<RANK>& frees, const Acts&... acts)
    : Base(frees,acts...), ExactBase(getFrees(),getActs()), HamiltonianBase(getFrees(),getActs()), LiouvilleanBase(getFrees(),getActs()) {}
    // Note that Frees and Acts are stored by value in Base, that’s why we need the getters here
 
};

#undef BASE_class


// The following provides a much more convenient interface:

namespace composite {

#define DISPATCHER(EX,HA,LI) (all(systemCharacteristics==SystemCharacteristics{EX,HA,LI})) return std::make_shared<Composite<EX,HA,LI,Acts...>>(frees,acts...)

template<typename... Acts>
std::shared_ptr<const composite::Base<Acts...>> make(const Acts&... acts)
{
  using namespace structure;
  
  const auto frees{[&]() {
    Frees<rank_v<Acts...>> res; res.fill({});
    
    hana::for_each(hana::make_tuple(acts...),[&](const auto& act) {
      int i=0;
      hana::for_each(act.retainedAxes, [&] (auto idx) {
        if (res[idx].get()) {
          if (res[idx].get()!=act->getFrees()[i]) throw composite::ConsistencyException(idx,i);
        }
        else res[idx]=SubSystemFree(act->getFrees()[i]);
        i++;
      });
    });
    
    return res;
  } ()};

  auto systemCharacteristics{hana::fold( hana::make_tuple(acts...), 
    hana::fold(tmptools::ordinals<Base<Acts...>::RANK>, SystemCharacteristics{false,false,false}, [&](SystemCharacteristics sc, auto idx) {
      const SubSystemFree& free=frees[idx];
      return sc || SystemCharacteristics{free.getEx()!=0,free.getHa()!=0,free.getLi()!=0};
    }), [&](SystemCharacteristics sc, const auto& act) {
      return sc || SystemCharacteristics{act.getEx()!=0,act.getHa()!=0,act.getLi()!=0};
    })};
  
  if      DISPATCHER(true ,true ,true ) ;
  else if DISPATCHER(true ,true ,false) ;
  else if DISPATCHER(true ,false,true ) ;
  else if DISPATCHER(true ,false,false) ;
  else if DISPATCHER(false,true ,true ) ;
  else if DISPATCHER(false,true ,false) ;
  else if DISPATCHER(false,false,true ) ;
  else return std::make_shared<Composite<false,false,false,Acts...>>(frees,acts...);
}


#undef DISPATCHER



} // composite



/// Template alias for backward compatibility
template<int... RetainedAxes>
using Act = composite::_<RetainedAxes...>;


#endif // CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED



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
 * \tparam IS_HA governs whether the class should inherit from composite::Hamiltonianhana::tuple<Acts...>
 * \tparam IS_LI governs whether the class should inherit from composite::Liouvillean
 */
