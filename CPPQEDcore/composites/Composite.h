// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED
#define   CPPQEDCORE_COMPOSITES_COMPOSITE_H_INCLUDED

#include "Interaction.h"
#include "Structure.h"

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
using Frees = typename ::structure::Interaction<RANK>::Frees;


/// Helper to create Composites
/**
 * It combines a set of \link retainedindexpositionsdefined retained axes\endlink (as tmptools::Vector) with a certain structure::Interaction element
 * of the corresponding arity (equalling the number of retained index positions) – this means that the given structure::Interaction element
 * will act on the given retained index positions.
 */
template<int... RetainedAxes>
struct _
{
public:
  explicit _(::structure::InteractionPtr<sizeof...(RetainedAxes)> ia) : ia_{ia} {}

  template <typename IA, typename ... T>
  static _ make(T&&... t)
  {
    return _{std::make_shared<const IA>(std::forward<T>(t)...)};
  }

  static constexpr auto retainedAxes=tmptools::vector<RetainedAxes...>;
  
  ::structure::InteractionPtr<sizeof...(RetainedAxes)> ia_;
  
};


template<typename... Acts>
constexpr auto rank_v=hana::maximum(hana::flatten(hana::make_tuple(Acts::retainedAxes...)))+1;



template<typename... Acts>
class Base
  : public ::structure::QuantumSystem<rank_v<Acts...>>,
    public ::structure::Averaged<rank_v<Acts...>>
{
public:
  // The calculated RANK
  static constexpr int RANK = rank_v<Acts...>;

protected:
  Base(Frees<RANK> frees, Acts... acts) : ::structure::QuantumSystem<RANK>([=] () { // calculate the dimensions
    typename structure::QuantumSystem<RANK>::Dimensions res;
    hana::for_each(tmptools::ordinals<RANK>,[=,&res](auto t) {res(t)=frees[t]->getDimension();});
    return res;
  }()), frees_{frees}, acts_{acts...} {}
  
  // Compile-time sanity check: each ordinal up to RANK has to be contained by at least one act.
  // (That is, no free system can be left out of the network of interactions)
  static_assert( hana::fold (tmptools::ordinals<RANK>, hana::true_c, [](auto state, auto ordinal) {
    return state && hana::contains(hana::flatten(hana::make_tuple(Acts::retainedAxes...)),ordinal);
  }) == true , "Composite not consistent" );

public:

  template<structure::LiouvillianAveragedTag LA>
  static std::ostream& streamKeyLA(std::ostream& os, size_t& i, Frees<RANK> frees, hana::tuple<Acts...> acts)
  {
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {::structure::streamKey<LA>(frees[idx],os,i);});
    hana::for_each(acts,[&](auto act) {::structure::streamKey<LA>(act.ia_,os,i);});    
    return os;
  }
  
  
  template<structure::LiouvillianAveragedTag LA>
  static size_t nAvrLA(Frees<RANK> frees, hana::tuple<Acts...> acts)
  {    
    return hana::fold(acts,
                      hana::fold(tmptools::ordinals<RANK>,
                                 0,
                                 [=](size_t res, auto idx) {return res + ::structure::nAvr<LA>(frees[idx]);}),
                      [=](size_t s, auto act) {return s+::structure::nAvr<LA>(act.ia_);});
  }
  

#ifndef NDEBUG
#pragma GCC warning "TODO: This solution is a bit insane, could be solved more effectively with blitz::Array slice"
#endif // NDEBUG
  template<structure::LiouvillianAveragedTag LA>
  static structure::Averages averageLA(double t, const quantumdata::LazyDensityOperator<RANK>& ldo, Frees<RANK> frees, hana::tuple<Acts...> acts, size_t numberAvr)
  {
    std::list<Averages> seqAverages{RANK+hana::size(acts).value}; // Averages res(numberAvr); size_t resIdx=0;
    
    {
      typename std::list<Averages>::iterator iter(seqAverages.begin());
      
      const auto lambda = [=,&iter,&ldo](auto av_or_li, auto v) {
        iter++->reference(quantumdata::partialTrace<std::decay_t<decltype(v)>>(ldo,[&](const auto& ldoS) {
          return structure::average<LA>(av_or_li,t,ldoS);
        }));
      };
      
      hana::for_each(tmptools::ordinals<RANK>,[=](auto idx) {lambda(frees[idx],tmptools::vector<idx>);});
      
      hana::for_each(acts,[=](auto act) {lambda(act.ia_,act.retainedAxes);});
      
    }
    
    Averages res(numberAvr); res=0;
    return cppqedutils::concatenate(seqAverages,res);
    
  }


private:
  // Implementing QuantumSystem interface

  double highestFrequency_v( ) const override
  {
    return std::max(
      frees_[hana::maximum(tmptools::ordinals<RANK>,[this](auto idx1, auto idx2) {
        return  frees_[idx1]->highestFrequency() < frees_[idx2]->highestFrequency();
      })]->highestFrequency(),
      hana::maximum( hana::transform( acts_, [](auto act) {return act.ia_->highestFrequency();} ) )
    );
  }
  
  std::ostream& streamParameters_v(std::ostream& os) const override
  {
    os<<"Composite\nDimensions: "<<this->getDimensions()<<". Total: "<<this->getTotalDimension()<<std::endl;
    
    hana::for_each(tmptools::ordinals<RANK>,[=,&os](auto idx) {
      os<<"Subsystem Nr. "<<idx<<std::endl;
      frees_[idx]->streamParameters(os);
    });
    
    hana::for_each(acts_, [=,&os](auto act) {
      hana::for_each(act.retainedAxes,[&](auto idx) {os<<idx<<" - ";});
      os<<"Interaction\n";
      act.ia_->streamParameters(os);
    });
    
    return os;
  }
  

  // Implementing the Averaged base

  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return streamKeyLA<::structure::LA_Av>(os,i, frees_ ,acts_ );}
  
  size_t nAvr_v() const override {return nAvrLA<::structure::LA_Av>( frees_,acts_ );}
  
  const ::structure::Averages average_v(double t, const quantumdata::LazyDensityOperator<RANK>& ldo) const override {return averageLA<::structure::LA_Av>(t,ldo,frees_,acts_,nAvr_v());}
  
  
  void process_v(::structure::Averages& avr) const override
  {
    ptrdiff_t l=-1, u=0;

    const auto lambda=[&](auto ptr) {
      if (const auto av=::structure::castAv(ptr); av && (u=l+av->nAvr())>l) {
        Averages temp(avr(blitz::Range(l+1,u)));
        av->process(temp);
        l=u;
      }
    };
  
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx]);});
    
    hana::for_each(acts_,[&](auto act) {lambda(act.ia_);});
  }
  
  
  std::ostream& stream_v(const structure::Averages& avr, std::ostream& os, int precision) const override
  {
    ptrdiff_t l=-1, u=0;
  
    const auto lambda=[&](auto ptr){
      if (const auto av=::structure::castAv(ptr); av && (u=l+av->nAvr())>l) {
        av->stream(avr(blitz::Range(l+1,u)),os,precision);
        l=u;
      }
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx]);});
    
    hana::for_each(acts_,[&](auto act) {lambda(act.ia_);});
    
    return os;
    
  }
  
  const Frees<RANK> frees_;
  
  const hana::tuple<Acts...> acts_;

};




template<typename... Acts>
class Exact
  : public structure::Exact<rank_v<Acts...>>
{
private:
  static constexpr int RANK=rank_v<Acts...>;

protected:
  Exact(Frees<RANK> frees, Acts... acts) : frees_{frees}, acts_{acts...} {}

private:
  void actWithU_v(double t, quantumdata::StateVectorLow<RANK>& psi, double t0) const override
  {
    const auto lambda=[&](auto ptr, auto v) {
      if (const auto ex=::structure::castEx(ptr)) for (auto& psiS : cppqedutils::sliceiterator::fullRange<decltype(v)>(psi)) ex->actWithU(t,psiS,t0);
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx],tmptools::vector<idx>);});
    
    hana::for_each(acts_,[&](auto act) {lambda(act.ia_,act.retainedAxes);});
    
  }
  
  
  bool applicableInMaster_v( ) const override
  {
    if (hana::fold(tmptools::ordinals<RANK>,true,[&](bool s, auto idx) {return s && ::structure::applicableInMaster(frees_[idx]);}))
      return hana::fold(acts_,true,[&](bool s, auto act) {return s && applicableInMaster(act.ia_);});
    else return false;
  }
  

  const Frees<RANK> frees_;
  const hana::tuple<Acts...> acts_;

};


template<typename... Acts>
class Hamiltonian
  : public structure::Hamiltonian<rank_v<Acts...>>
{
private:
  static constexpr int RANK=rank_v<Acts...>;

protected:
  Hamiltonian(Frees<RANK> frees, Acts... acts) : frees_{frees}, acts_{acts...} {}

private:
  void addContribution_v(double t, const quantumdata::StateVectorLow<RANK>& psi, quantumdata::StateVectorLow<RANK>& dpsidt, double t0) const override
  {
    const auto lambda=[&](auto ptr, auto v) {
      if (const auto ha=::structure::castHa(ptr))
        boost::for_each(cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(psi),
                        cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(dpsidt),
                        [&](const auto& psiS, auto& dpsidtS) {ha->addContribution(t,psiS,dpsidtS,t0);});
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx],tmptools::vector<idx>);});

    hana::for_each(acts_,[&](auto act) {lambda(act.ia_,act.retainedAxes);});
    
  }
  
  const Frees<RANK> frees_;
  const hana::tuple<Acts...> acts_;

};


template<typename... Acts>
class Liouvillian
  : public structure::Liouvillian<rank_v<Acts...>>
{
private:
  static constexpr int RANK=rank_v<Acts...>;

protected:
  Liouvillian(Frees<RANK> frees, Acts... acts) : frees_{frees}, acts_{acts...} {}

private:
  std::ostream& streamKey_v(std::ostream& os, size_t& i) const override {return Base<Acts...>::template streamKeyLA<structure::LA_Li>(os,i, frees_,acts_);}
  
  size_t nAvr_v() const override {return Base<Acts...>::template nAvrLA<structure::LA_Li>(frees_,acts_);}
  
  const Rates average_v(double t, const quantumdata::LazyDensityOperator<RANK>& ldo) const override
  {
    return Base<Acts...>::template averageLA<structure::LA_Li>(t,ldo,frees_,acts_,nAvr_v());
  }

  void actWithJ_v(double t, quantumdata::StateVectorLow<RANK>& psi, size_t ordoJump) const override
  {
    bool flag=false;
    
    const auto lambda=[&](auto ptr, auto v) {
      if (const auto li=::structure::castLi(ptr); !flag && li) {
        size_t n=li->nAvr();
        if (ordoJump<n) {
          for (auto& psiS : cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(psi)) li->actWithJ(t,psiS,ordoJump);
          flag=true;
        }
        ordoJump-=n;  
      }
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx],tmptools::vector<idx>);});
    
    hana::for_each(acts_,[&](const auto& act) {lambda(act.ia_,act.retainedAxes);});
    
  }
  
  
  void actWithSuperoperator_v(double t, const quantumdata::DensityOperatorLow<RANK>& rho, quantumdata::DensityOperatorLow<RANK>& drhodt, size_t ordoJump) const override
  {
    bool flag=false;
    
    const auto lambda=[&](auto ptr, auto ev) {
      if (const auto li=::structure::castLi(ptr); !flag && li) {
        size_t n=li->nAvr();
        if (ordoJump<n) {
          using ExtendedV=std::decay_t<decltype(ev)>;

          boost::for_each(cppqedutils::sliceiterator::fullRange<ExtendedV>(rho),cppqedutils::sliceiterator::fullRange<ExtendedV>(drhodt),
                          [&](const auto& rhoS, auto& drhodtS) {li->actWithSuperoperator(t,rhoS,drhodtS,ordoJump);});
          
          flag=true;
        }
        ordoJump-=n;  
      }
    };
    
    hana::for_each(tmptools::ordinals<RANK>,[&](auto idx) {lambda(frees_[idx],tmptools::vector<idx,idx+RANK>);});
    
    hana::for_each(acts_,[&](auto act) {lambda(act.ia_,hana::concat(act.retainedAxes,
                                                                    hana::transform(act.retainedAxes,
                                                                                    [](auto element) constexpr {return element+hana::int_c<RANK>;} ) ) ); } );
    
  }
  
  const Frees<RANK> frees_;
  const hana::tuple<Acts...> acts_;

};


template<typename TAG>
class EmptyBase
{
public:
  template<typename... T>
  EmptyBase(T...) {}
  
};


} // composite



#define BASE_class(Aux,Class) std::conditional_t<IS_##Aux,composite::Class<Acts...>,composite::EmptyBase<composite::Class<Acts...>>>


template<bool IS_EX, bool IS_HA, bool IS_LI, typename... Acts>
// Acts should model a fusion sequence of Acts
class Composite
  : public composite::Base<Acts...>,
    public BASE_class(EX,Exact),
    public BASE_class(HA,Hamiltonian),
    public BASE_class(LI,Liouvillian)
{
public:
  typedef composite::Base<Acts...> Base;
  
  typedef BASE_class(EX,Exact) ExactBase;
  typedef BASE_class(HA,Hamiltonian) HamiltonianBase;
  typedef BASE_class(LI,Liouvillian) LiouvillianBase;
  
  // Constructor
  Composite(composite::Frees<Base::RANK> frees, Acts... acts)
    : Base(frees,acts...), ExactBase(frees,acts...), HamiltonianBase(frees,acts...), LiouvillianBase(frees,acts...) {}
    // Note that Frees and Acts are stored by value in Base, that’s why we need the getters here
 
};

#undef BASE_class


// The following provides a much more convenient interface:

namespace composite {

#define DISPATCHER(EX,HA,LI) (all(sysChar==SystemCharacteristics{EX,HA,LI})) return std::make_shared<Composite<EX,HA,LI,Acts...>>(frees,acts...)

template<typename... Acts>
std::shared_ptr<const composite::Base<Acts...>> make(Acts... acts)
{
  using namespace structure;
  
  const auto frees{[&]() {
    Frees<rank_v<Acts...>> res;// res.fill({});
    
    hana::for_each(hana::make_tuple(acts...),[&](auto act) {
      int i=0;
      hana::for_each(act.retainedAxes, [&] (auto idx) {
        if (res[idx]) {
          if ( res[idx]!=(*act.ia_)[i] ) throw composite::ConsistencyException(idx,i);
        }
        else res[idx]=(*act.ia_)[i];
        i++;
      });
    });
    
    return res;
  } ()};
  
  SystemCharacteristics sysChar{hana::fold( hana::make_tuple(acts...), 
    hana::fold(tmptools::ordinals<rank_v<Acts...>>, SystemCharacteristics{false,false,false}, [&](SystemCharacteristics sc, auto idx) {
      const auto freePtr=frees[idx];
      return sc or systemCharacteristics(freePtr);
    }), [&](SystemCharacteristics sc, auto act) {
      return sc or systemCharacteristics(act.ia_);
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
 * \tparam IS_LI governs whether the class should inherit from composite::Liouvillian
 */
