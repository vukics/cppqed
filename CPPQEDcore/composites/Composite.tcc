// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_COMPOSITES_COMPOSITE_TCC_INCLUDED
#define   CPPQEDCORE_COMPOSITES_COMPOSITE_TCC_INCLUDED

#include "Composite.h"

#include "Interaction.h"
#include "LazyDensityOperator.h"

#include "Algorithm.h"
#include "SliceIterator.tcc"

#include <boost/mpl/for_each.hpp>

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/algorithm/iteration/fold.hpp>

#include <boost/range/algorithm/max_element.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include <algorithm>
#include <list>



// template<int>
struct CompositeConsistencyException : std::logic_error
{
  CompositeConsistencyException(int idxx, int ii) : std::logic_error(std::to_string(idxx)+" "+std::to_string(ii)), idx(idxx), i(ii) {}

  const int idx, i;
};


///////////////
//
// Construction
//
///////////////


template<typename VA>
auto composite::fillFrees(const VA& acts)
{
  typename composite::Base<VA>::Frees res; res.fill(SubSystemFree());
  
  boost::fusion::for_each(acts,[&](const auto& act) -> void {
    int i=0;
    mpl::for_each<std::decay_t<decltype(act)>>([&](auto t) -> void {
      static const int idx=decltype(t)::value;
      if (res[idx].get()) {
        if (res[idx].get()!=act.get()->getFrees()[i]) throw CompositeConsistencyException(idx,i);
      }
      else res[idx]=SubSystemFree(act.get()->getFrees()[i]);
      i++;
    });
  });

  return res;
}


template<int RANK>
auto composite::RankedBase<RANK>::fillDimensions(const Frees& frees)
{
  typename composite::RankedBase<RANK>::Dimensions res;
  
  mpl::for_each<Ordinals>([&](auto t) -> void {
    static const int idx=decltype(t)::value;
    res(idx)=frees[idx].get()->getDimension();
  });
  
  return res;
}



/////////////
//
// Parameters
//
/////////////

template<typename VA>
std::ostream& composite::Base<VA>::streamParameters_v(std::ostream& os) const
{
  os<<"Composite\nDimensions: "<<RBase::getDimensions()<<". Total: "<<RBase::getTotalDimension()<<std::endl;
  
  mpl::for_each<Ordinals>([&](auto t)-> void {
    static const int idx=decltype(t)::value;
    os<<"Subsystem Nr. "<<idx<<std::endl;
    frees_[idx].get()->streamParameters(os);
  });
  
  boost::fusion::for_each(
    acts_,
    [&](const auto& act) -> void {
      using Act=std::decay_t<decltype(act)>;
      mpl::for_each<Act>([&](auto t) {
        os<<decltype(t)::value<<" - ";
      });
      os<<"Interaction\n";
      act.get()->streamParameters(os);
    });
  
  return os;
}

//////////////
//
// Frequencies
//
//////////////


template<typename VA>
double composite::Base<VA>::highestFrequency_v() const
{
  return boost::max_element(
    frees_,
    [](const SubSystemFree& ssf1, const SubSystemFree& ssf2) {
      return ssf1.get()->highestFrequency() < ssf2.get()->highestFrequency();
    })->get()->highestFrequency();
  
  // NEEDS_WORK add the interactions here
}


////////
//
// Exact
//
////////


template<typename VA>
bool composite::Exact<VA>::applicableInMaster_v() const
{
  bool res=true;
  {
    mpl::for_each<Ordinals>([&](auto t) -> void {
      res&=frees_[decltype(t)::value].applicableInMaster();
    });
    
    if (res)
      res=boost::fusion::fold(acts_,res,[&](bool s, const auto& act) -> bool {
        return s && act.applicableInMaster();
      });
  }
  return res;
}



template<typename VA>
void composite::Exact<VA>::actWithU_v(double t, StateVectorLow& psi, double t0) const
{
  const auto lambda=[&](auto ex, auto v) -> void {if (ex) for (auto& psiS : cppqedutils::sliceiterator::fullRange<decltype(v)>(psi)) ex->actWithU(t,psiS,t0);};

  mpl::for_each<Ordinals>([&](auto t) -> void {
    static const int idx=decltype(t)::value;
    lambda(frees_[idx].getEx(),tmptools::Vector<idx>());
  });
  boost::fusion::for_each(acts_,[&](const auto& act) -> void {
    lambda(act.getEx(),act);
  });

}


//////////////
//
// Hamiltonian
//
//////////////

template<typename VA>
void composite::Hamiltonian<VA>::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const
{
  const auto lambda=[&](auto ha, auto v) -> void {
    using Vec=std::decay_t<decltype(v)>;
    if (ha)
      boost::for_each(cppqedutils::sliceiterator::fullRange<Vec>(psi),
                      cppqedutils::sliceiterator::fullRange<Vec>(dpsidt),
                      [&](const auto& psiS, auto& dpsidtS) {ha->addContribution(t,psiS,dpsidtS,t0);}); 
  };
  
  mpl::for_each<Ordinals>([&](auto t) -> void {
    static const int idx=decltype(t)::value;
    lambda(frees_[idx].getHa(),tmptools::Vector<idx>());
  });
  boost::fusion::for_each(acts_,[&](const auto& act) -> void {
    lambda(act.getHa(),act);
  });
  
}


///////////////////////
//
// Liouvillian-Averaged
//
///////////////////////

template<typename VA> template<structure::LiouvillianAveragedTag LA>
std::ostream& composite::Base<VA>::streamKeyLA(std::ostream& os, size_t& i, const Frees& frees, const VA& acts)
{
  mpl::for_each<Ordinals>([&](auto t) -> void {
    static const int idx=decltype(t)::value;
    frees[idx].template streamKey<LA>(os,i);
  });
  boost::fusion::for_each(acts,[&](const auto& act) -> void {
    act.template streamKey<LA>(os,i);
  });
  
  return os;
}


template<typename VA> template<structure::LiouvillianAveragedTag LA>
size_t composite::Base<VA>::nAvrLA(const Frees& frees, const VA& acts)
{
  size_t res=0;

  mpl::for_each<Ordinals>([&](auto t) -> void {
    res+=frees[decltype(t)::value].template nAvr<LA>();
  });

  return boost::fusion::fold(acts,res,[&](size_t s, const auto& act) -> size_t {
    return s+act.template nAvr<LA>();
  });
  
}



template<typename VA> template<structure::LiouvillianAveragedTag LA>
auto
composite::Base<VA>::averageLA(double t, const LazyDensityOperator& ldo, const Frees& frees, const VA& acts, size_t numberAvr) -> const Averages
{
  std::list<Averages> seqAverages(RANK+mpl::size<VA>::value);

  {
    typename std::list<Averages>::iterator iter(seqAverages.begin());

    const auto lambda=[&](auto av, auto v) {
      iter++->reference(quantumdata::partialTrace<std::decay_t<decltype(v)>>(ldo,[&](const auto& ldoS) {
        return structure::average(av,t,ldoS);
      }));
    };

    auto tag=structure::LiouvillianAveragedTag_<LA>();
    
    mpl::for_each<Ordinals>([&](auto t) -> void {
      static const int idx=decltype(t)::value;
      lambda(frees[idx].getLA(tag),tmptools::Vector<idx>());
    });
    
    boost::fusion::for_each(acts,[&](const auto& act) -> void {
      lambda(act.getLA(tag),act);
    });
    
  }

  Averages res(numberAvr); res=0;
  return cppqedutils::concatenate(seqAverages,res);

}


//////////////
//
// Liouvillian
//
//////////////


template<typename VA>
void composite::Liouvillian<VA>::actWithJ_v(double t, StateVectorLow& psi, size_t ordoJump) const
{
  bool flag=false;

  const auto lambda=[&](auto li, auto v) -> void {
    if (!flag && li) {
      size_t n=li->nAvr();
      if (ordoJump<n) {
        for (auto& psiS : cppqedutils::sliceiterator::fullRange<std::decay_t<decltype(v)>>(psi)) li->actWithJ(t,psiS,ordoJump);
        flag=true;
      }
      ordoJump-=n;  
    }
  };

  mpl::for_each<Ordinals>([&](auto t) -> void {
    static const int idx=decltype(t)::value;
    lambda(frees_[idx].getLi(),tmptools::Vector<idx>());
  });
  boost::fusion::for_each(acts_,[&](const auto& act) -> void {
    lambda(act.getLi(),act);
  });

}


template<typename VA>
void composite::Liouvillian<VA>::actWithSuperoperator_v(double t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t ordoJump) const
{
  bool flag=false;

  const auto lambda=[&](auto li, auto v) -> void {
    using Vec=std::decay_t<decltype(v)>;
    if (!flag && li) {
      size_t n=li->nAvr();
      if (ordoJump<n) {
        boost::for_each(cppqedutils::sliceiterator::fullRange<tmptools::ExtendVector_t<RANK,Vec>>(rho),
                        cppqedutils::sliceiterator::fullRange<tmptools::ExtendVector_t<RANK,Vec>>(drhodt),
                        [&](const auto& rhoS, auto& drhodtS) {li->actWithSuperoperator(t,rhoS,drhodtS,ordoJump);});
        flag=true;
      }
      ordoJump-=n;  
    }
  };

  mpl::for_each<Ordinals>([&](auto t) -> void {
    static const int idx=decltype(t)::value;
    lambda(frees_[idx].getLi(),tmptools::Vector<idx>());
  });
  boost::fusion::for_each(acts_,[&](const auto& act) -> void {
    lambda(act.getLi(),act);
  });

}


///////////
//
// Averaged
//
///////////


template<typename VA>
void
composite::Base<VA>::process_v(Averages& avr) const
{
  ptrdiff_t l=-1, u=0;

  const auto lambda=[&](auto av) -> void {
    using blitz::Range;
    if (av && (u=l+av->nAvr())>l) {
      Averages temp(avr(Range(l+1,u)));
      av->process(temp);
      l=u;
    }
  };

  mpl::for_each<Ordinals>([&](auto t) -> void {lambda(frees_[decltype(t)::value].getAv());});
  boost::fusion::for_each(acts_,[&](const auto& act) -> void {lambda(act.getAv());});

}



template<typename VA>
std::ostream&
composite::Base<VA>::stream_v(const Averages& avr, std::ostream& os, int precision) const
{
  ptrdiff_t l=-1, u=0;
  
  const auto lambda=[&](auto av) -> void {
    using blitz::Range;
    if (av && (u=l+av->nAvr())>l) {
      av->stream(avr(Range(l+1,u)),os,precision);
      l=u;
    }
  };
  
  mpl::for_each<Ordinals>([&](auto t) -> void {lambda(frees_[decltype(t)::value].getAv());});
  boost::fusion::for_each(acts_,[&](const auto& act) -> void {lambda(act.getAv());});

  return os;
  
}



////////
//
// Maker
//
////////


#define DISPATCHER(EX,HA,LI) (all(systemCharacteristics==SystemCharacteristics{EX,HA,LI})) return std::make_shared<Composite<VA,EX,HA,LI> >(frees,acts)

template<typename VA>
const typename composite::Base<VA>::Ptr composite::doMake(const VA& acts)
{
  using namespace structure;
  
  const typename Base<VA>::Frees frees(fillFrees(acts));

  SystemCharacteristics systemCharacteristics{false,false,false};
  
  mpl::for_each<typename Base<VA>::Ordinals>([&](auto t) -> void {
    const SubSystemFree& free=frees[std::decay_t<decltype(t)>::value];
    systemCharacteristics|=SystemCharacteristics{free.getEx()!=0,free.getHa()!=0,free.getLi()!=0};
  });
  
  systemCharacteristics=boost::fusion::fold(acts,systemCharacteristics,[&](SystemCharacteristics sc, const auto& act) -> SystemCharacteristics {
    return sc || SystemCharacteristics{act.getEx()!=0,act.getHa()!=0,act.getLi()!=0};
  });
    
  if      DISPATCHER(true ,true ,true ) ;
  else if DISPATCHER(true ,true ,false) ;
  else if DISPATCHER(true ,false,true ) ;
  else if DISPATCHER(true ,false,false) ;
  else if DISPATCHER(false,true ,true ) ;
  else if DISPATCHER(false,true ,false) ;
  else if DISPATCHER(false,false,true ) ;
  else return std::make_shared<Composite<VA,false,false,false> >(frees,acts);
}


#undef DISPATCHER

#undef CALL_composite_worker
#endif // CPPQEDCORE_COMPOSITES_COMPOSITE_TCC_INCLUDED
