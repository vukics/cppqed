// -*- C++ -*-
#ifndef COMPOSITE_AVERAGED_LIOUVILLEAN_ITERATING

#ifndef   ELEMENTS_COMPOSITES_IMPL_COMPOSITE_TCC_INCLUDED
#define   ELEMENTS_COMPOSITES_IMPL_COMPOSITE_TCC_INCLUDED

#include "Composite.h"

#include "Interaction.h"

#include "impl/LazyDensityOperator.tcc"

#include "Exception.h"

#include "Algorithm.h"
#include "Range.h"

#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/algorithm/iteration/fold.hpp>

#include <boost/mpl/for_each.hpp>

#include <boost/bind.hpp>

#include <boost/lambda/lambda.hpp>
namespace bll=boost::lambda;

#include <algorithm>
#include <list>


#define CALL_composite_worker(object) composite::worker<Ordinals>(acts_,object);


// NEEDS_WORK in the helper classes worker member functions could be factored out. (Eg ActWithU --- the amount of code we save is incredible!)
// Similarly, like in BinarySystem, there is too much boilerplate.


// template<int>
struct CompositeConsistencyException : cpputils::Exception
{
  CompositeConsistencyException(int idxx, int ii) : idx(idxx), i(ii) {}

  const int idx, i;
};


//////////////////
//
// Overall helpers
//
//////////////////


namespace composite {


template<typename Ordinals, typename VA, typename H>
void worker(const VA& acts, const H& helper)
{
  mpl::for_each<Ordinals>(helper);
  boost::fusion::for_each(acts,helper);

  // IMPORTANT NOTICE!!! Both mpl::for_each and boost::fusion::for_each stores the helper BY VALUE. So, what happens is that the initial helper is COPIED to both, so the changes made during mpl::for_each in those members of the helper which are stored BY VALUE in H, will be "undone" at startup of boost::fusion::for_each. Ergo, there MUST NOT BE such members. (There can be const members stored by value though, of course.) ... How about static members then?

}


///////////////
//
// Construction
//
///////////////



bool compareFreesFrequency(const SubSystemFree& ssf1, const SubSystemFree& ssf2);


template<int RANK>
class FillFrees
{
public:
  typedef blitz::TinyVector<SubSystemFree,RANK> Frees;

  FillFrees(Frees& frees) : frees_(frees) {}

  
  template<typename Act>
  class Inner
  {
  public:
    Inner(const Act& act, Frees& frees) : act_(act), frees_(frees), i_(0) {}

    template<typename T>
    void operator()(T)
    {
      static const int idx=T::value;
      if (frees_(idx).get()) {
	if (frees_(idx).get()!=act_.get()->getFrees()(i_)) throw CompositeConsistencyException(idx,i_);
      }
      else frees_(idx)=SubSystemFree(act_.get()->getFrees()(i_));
      i_++;
    }

  private:
    const Act& act_;
    Frees& frees_;
    
    int i_;
    
  };


  template<typename Act>
  void operator()(const Act& act) const
  {
    mpl::for_each<Act>(Inner<Act>(act,frees_));
  }

private:
  Frees& frees_;

};


template<int RANK>
class FillDimensions
{
public:
  typedef blitz::TinyVector<SubSystemFree,RANK> Frees;
  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;

  FillDimensions(const Frees& frees, Dimensions& dims) : frees_(frees), dims_(dims) {}

  template<typename T>
  void operator()(T) const
  {
    dims_(T::value)=frees_(T::value).get()->getDimension();
  }

private:
  const Frees& frees_;
  Dimensions&   dims_;

};


} // composite

#define RETURN_type typename composite::Base<VA>::Frees

template<typename VA>
const RETURN_type
composite::fillFrees(const VA& acts)
{
  RETURN_type res; res=SubSystemFree();
  boost::fusion::for_each(acts,FillFrees<blitzplusplus::TinyVectorLengthTraits<RETURN_type>::value>(res));
  return res;
}

#undef  RETURN_type
#define RETURN_type typename composite::RankedBase<RANK>::Dimensions

template<int RANK>
const RETURN_type
composite::RankedBase<RANK>::fillDimensions(const Frees& frees)
{
  RETURN_type res;
  mpl::for_each<Ordinals>(FillDimensions<RANK>(frees,res));
  return res;
}

#undef  RETURN_type



/////////////
//
// Parameters
//
/////////////



template<typename VA>
class composite::Base<VA>::DisplayParameters
{
public:
  DisplayParameters(const Frees& frees, std::ostream& os) : frees_(frees), os_(os) {}


  template<typename Act>
  class Inner
  {
  public:
    Inner(const Act& act, std::ostream& os) : act_(act), os_(os) {}

    template<typename T>
    void operator()(T)
    {
      os_<<T::value<<" - ";
    }

  private:
    const Act& act_;
    std::ostream& os_;
    
  };


  template<typename Act>
  void operator()(const Act& act) const
  {
    os_<<"# ";
    mpl::for_each<Act>(Inner<Act>(act,os_));
    os_<<"Interaction\n";
    act.get()->displayParameters(os_);
  }


  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    os_<<"# Subsystem Nr. "<<IDX<<std::endl;
    frees_(IDX).get()->displayParameters(os_);
  }
  
private:
  const Frees& frees_;

  std::ostream& os_;

};



template<typename VA>
void composite::Base<VA>::displayParameters_v(std::ostream& os) const
{
  os<<"# Composite\n# Dimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<std::endl;
  CALL_composite_worker( DisplayParameters(frees_,os) ) ;
}

//////////////
//
// Frequencies
//
//////////////


template<typename VA>
double composite::Base<VA>::highestFrequency_v() const
{
  return boost::max_element(frees_,composite::compareFreesFrequency)->get()->highestFrequency();
  // NEEDS_WORK add the interactions here
}


////////
//
// Exact
//
////////


template<typename VA>
class composite::Exact<VA>::IsUnitary
{
public:
  IsUnitary(const Frees& frees, bool& isIt) : frees_(frees), isIt_(isIt) {}

  typedef bool result_type;
  
  template<typename Act>
  bool operator()(bool s, const Act& act)
  {
    return s && act.isUnitary();
  }

  template<typename T>
  void operator()(T) const
  {
    isIt_&=frees_(T::value).isUnitary();
  }

private:
  const Frees& frees_;

  bool& isIt_;

};



template<typename VA>
bool composite::Exact<VA>::isUnitary_v() const
{
  bool res=true;
  {
    const IsUnitary helper(frees_,res);
    mpl::for_each<Ordinals>(helper);
    if (res) res=boost::fusion::fold(acts_,res,helper);
  }
  return res;
}



template<typename VA>
class composite::Exact<VA>::ActWithU
{
public:
  ActWithU(const Frees& frees, double dtdid, StateVectorLow& psi) : frees_(frees), dtdid_(dtdid), psi_(psi) {}

  template<typename Vec, typename Ex>
  void help(typename Ex::Ptr ex) const
  {
    if (ex) boost::for_each(blitzplusplus::basi::fullRange<Vec>(psi_),boost::bind(&Ex::actWithU,ex,dtdid_,::_1));
    // namespace qualification :: is necessary because otherwise :: and mpl:: are equally good matches.
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<Act,typename Act::Ex>(act.getEx());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<tmptools::Vector<IDX>,composite::SubSystemFree::Ex>(frees_(IDX).getEx());
  }


private:
  const Frees& frees_;

  const double dtdid_;
  StateVectorLow& psi_; 

};


template<typename VA>
void composite::Exact<VA>::actWithU_v(double dtdid, StateVectorLow& psi) const
{
  CALL_composite_worker( ActWithU(frees_,dtdid,psi) ) ;
}


//////////////
//
// Hamiltonian
//
//////////////

template<typename VA>
class composite::Hamiltonian<VA>::AddContribution
{
public:
  AddContribution(const Frees& frees, double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) 
    : frees_(frees), t_(t), psi_(psi), dpsidt_(dpsidt), tIntPic0_(tIntPic0) {}

  template<typename Vec, typename Ha>
  void help(typename Ha::Ptr ha) const
  {
    if (ha) 
      cpputils::for_each(blitzplusplus::basi::fullRange<Vec>(psi_),
			 blitzplusplus::basi::begin<Vec>(dpsidt_),
			 boost::bind(&Ha::addContribution,ha,t_,::_1,::_2,tIntPic0_)); 
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<Act,typename Act::Ha>(act.getHa());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<tmptools::Vector<IDX>,composite::SubSystemFree::Ha>(frees_(IDX).getHa());
  }

private:
  const Frees& frees_;

  const double t_;
  const StateVectorLow& psi_; 
  StateVectorLow& dpsidt_; 
  const double tIntPic0_;

};


template<typename VA>
void composite::Hamiltonian<VA>::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  CALL_composite_worker( AddContribution(frees_,t,psi,dpsidt,tIntPic0) ) ;
}


//////////////
//
// Liouvillean
//
//////////////


template<typename VA>
class composite::Liouvillean<VA>::Probas
{
public:
  typedef typename std::list<Probabilities>::iterator Iter;

  Probas(const Frees& frees, double t, const LazyDensityOperator& ldo, Iter& iter) 
    : frees_(frees), t_(t), ldo_(ldo), iter_(iter) {}

  template<typename Vec, typename Li>
  void help(typename Li::Ptr li) const
  {
    iter_++->reference(quantumdata::partialTrace<Vec,Probabilities>(ldo_,boost::bind(structure::probabilities<Li::N_RANK>,li,t_,::_1)));    
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<Act,typename Act::Li>(act.getLi());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<tmptools::Vector<IDX>,composite::SubSystemFree::Li>(frees_(IDX).getLi());
  }

private:
  const Frees& frees_;

  const double t_;
  const LazyDensityOperator& ldo_;
  
  Iter& iter_;

};


template<typename VA>
const typename composite::Liouvillean<VA>::Probabilities
composite::Liouvillean<VA>::probabilities_v(double t, const LazyDensityOperator& ldo) const
{
  std::list<Probabilities> seqProbabilities(RANK+mpl::size<VA>::value);

  {
    typename Probas::Iter iter(seqProbabilities.begin());

    CALL_composite_worker( Probas(frees_,t,ldo,iter) ) ;
  }

  Probabilities res(nJumps_v()); res=0;
  return cpputils::concatenate(seqProbabilities,res);

}



template<typename VA>
class composite::Liouvillean<VA>::ActWithJ
{
public:
  ActWithJ(const Frees& frees, double t, StateVectorLow& psi, size_t& ordoJump, bool& flag) : frees_(frees), t_(t), psi_(psi), ordoJump_(ordoJump), flag_(flag) {}

  template<typename Vec, typename Li>
  void help(typename Li::Ptr li) const
  {
    if (!flag_ && li) {
      size_t n=li->nJumps();
      if (ordoJump_<n) {
	boost::for_each(blitzplusplus::basi::fullRange<Vec>(psi_),boost::bind(&Li::actWithJ,li,t_,::_1,ordoJump_));
	flag_=true;
      }
      ordoJump_-=n;  
    }
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<Act,typename Act::Li>(act.getLi());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<tmptools::Vector<IDX>,composite::SubSystemFree::Li>(frees_(IDX).getLi());
  }

private:
  const Frees& frees_;

  const double t_;
  StateVectorLow& psi_; 
  size_t& ordoJump_;
  bool& flag_;

};


template<typename VA>
void composite::Liouvillean<VA>::actWithJ_v(double t, StateVectorLow& psi, size_t ordoJump) const
{
  bool flag=false;
  CALL_composite_worker( ActWithJ(frees_,t,psi,ordoJump,flag) ) ;
}


#define COMPOSITE_AVERAGED_LIOUVILLEAN_ITERATING

#define A_or_L Liouvillean
#define Class Liouvillean
#define NumberName Jumps

#include "Composite.tcc"


  

///////////
//
// Averaged
//
///////////


#define COMPOSITE_AVERAGED_LIOUVILLEAN_ITERATING

#define A_or_L Averaged
#define Class  Base
#define NumberName Avr

#include "Composite.tcc"



template<typename VA>
class composite::Base<VA>::Average
{
public:

  typedef typename std::list<Averages>::iterator Iter;

  Average(const Frees& frees, double t, const LazyDensityOperator& ldo, Iter& iter) 
    : frees_(frees), t_(t), ldo_(ldo), iter_(iter) {}

  template<typename Vec, typename Av>
  void help(typename Av::Ptr av) const
  {
    iter_++->reference(quantumdata::partialTrace<Vec,Averages>(ldo_,boost::bind(structure::average<Av::N_RANK>,av,t_,::_1)));
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<Act,typename Act::Av>(act.getAv());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<tmptools::Vector<IDX>,composite::SubSystemFree::Av>(frees_(IDX).getAv());
  }

private:
  const Frees& frees_;

  const double t_;
  const LazyDensityOperator& ldo_;
  
  Iter& iter_;

};


template<typename VA>
const typename composite::Base<VA>::Averages
composite::Base<VA>::average_v(double t, const LazyDensityOperator& ldo) const
{
  std::list<Averages> seqAverages(RANK+mpl::size<VA>::value);

  {
    typename std::list<Averages>::iterator iter(seqAverages.begin());

    CALL_composite_worker( Average(frees_,t,ldo,iter) ) ;
  }

  Averages res(nAvr_v()); res=0;
  return cpputils::concatenate(seqAverages,res);

}




template<typename VA>
class composite::Base<VA>::Process
{
public:
  Process(const Frees& frees, const Averages& avr, ptrdiff_t& l, ptrdiff_t& u) 
    : frees_(frees), avr_(avr), l_(l), u_(u) {}


  template<typename Av>
  void help(typename Av::Ptr av) const
  {
    using blitz::Range;
    if (av && (u_=l_+av->nAvr())>l_) {
      Averages temp(avr_(Range(l_+1,u_)));
      av->process(temp);
    }
    std::swap(l_,u_);
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<typename Act::Av>(act.getAv());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<composite::SubSystemFree::Av>(frees_(IDX).getAv());
  }

private:
  const Frees& frees_;

  const Averages& avr_;

  ptrdiff_t &l_, &u_;

};


template<typename VA>
void
composite::Base<VA>::process_v(Averages& avr) const
{
  ptrdiff_t l=-1, u=0;
  CALL_composite_worker( Process(frees_,avr,l,u) ) ;
}



template<typename VA>
class composite::Base<VA>::Display
{
public:
  Display(const Frees& frees, const Averages& avr, std::ostream& os, int precision, ptrdiff_t& l, ptrdiff_t& u) 
    : frees_(frees), avr_(avr), os_(os), precision_(precision), l_(l), u_(u) {}

  template<typename Av>
  void help(typename Av::Ptr av) const
  {
    using blitz::Range;
    if (av && (u_=l_+av->nAvr())>l_) av->display(avr_(Range(l_+1,u_)),os_,precision_);
    std::swap(l_,u_);
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help<typename Act::Av>(act.getAv());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help<composite::SubSystemFree::Av>(frees_(IDX).getAv());
  }

private:
  const Frees& frees_;

  const Averages& avr_;
  std::ostream& os_;

  const int precision_;

  ptrdiff_t &l_, &u_;

};


template<typename VA>
void
composite::Base<VA>::display_v(const Averages& avr, std::ostream& os, int precision) const
{
  ptrdiff_t l=-1, u=0;
  CALL_composite_worker( Display(frees_,avr,os,precision,l,u) ) ;
}



////////
//
// Maker
//
////////


namespace composite {

template<typename VA>
class QuerySystemCharacteristics
{
public:
  typedef typename Base<VA>::Frees Frees;
  
  typedef structure::SystemCharacteristics result_type;

  QuerySystemCharacteristics(const Frees& frees, result_type& sc) : frees_(frees), sc_(sc) {}

  
  template<typename Act>
  result_type operator()(result_type sc, const Act& act)
  {
    const typename Act::QuantumSystemPtr ia=act.getQS();
    return sc || result_type(qse(ia),qsh(ia),qsl(ia));
  }

  template<typename T>
  void operator()(T) const
  {
    const structure::QuantumSystem<1>::Ptr free=frees_(T::value).get();
    sc_|=result_type(qse(free),qsh(free),qsl(free));
  }

private:
  const Frees& frees_;

  result_type& sc_;

};

} // composite


#define DISPATCHER(EX,HA,LI) (all(systemCharacteristics==SystemCharacteristics(EX,HA,LI))) return boost::make_shared<Composite<VA,EX,HA,LI> >(frees,acts)

template<typename VA>
const typename composite::Base<VA>::Ptr composite::doMake(const VA& acts)
{
  using namespace structure;
  
  const typename Base<VA>::Frees frees(fillFrees(acts));

  SystemCharacteristics systemCharacteristics(false,false,false);
  
  {
    const composite::QuerySystemCharacteristics<VA> helper(frees,systemCharacteristics);
    mpl::for_each<typename Base<VA>::Ordinals>(helper);
    systemCharacteristics=boost::fusion::fold(acts,systemCharacteristics,helper);
  }
    
  if      DISPATCHER(true ,true ,true ) ;
  else if DISPATCHER(true ,true ,false) ;
  else if DISPATCHER(true ,false,true ) ;
  else if DISPATCHER(true ,false,false) ;
  else if DISPATCHER(false,true ,true ) ;
  else if DISPATCHER(false,true ,false) ;
  else if DISPATCHER(false,false,true ) ;
  else return boost::make_shared<Composite<VA,false,false,false> >(frees,acts);
}


#undef DISPATCHER






#undef CALL_composite_worker

#endif // ELEMENTS_COMPOSITES_IMPL_COMPOSITE_TCC_INCLUDED


#else  // COMPOSITE_AVERAGED_LIOUVILLEAN_ITERATING


#define KEY_name BOOST_PP_CAT(A_or_L,Key)

namespace composite {

template<typename VA>
class BOOST_PP_CAT(Display,KEY_name)
{
public:
  typedef typename Base<VA>::Frees Frees;
  
  BOOST_PP_CAT(Display,KEY_name)(const Frees& frees, std::ostream& os, size_t& i) 
    : frees_(frees), os_(os), i_(i) {}

  template<typename Act>
  void operator()(const Act& act) const
  {
    act.BOOST_PP_CAT(display,KEY_name)(os_,i_);
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    frees_(IDX).BOOST_PP_CAT(display,KEY_name)(os_,i_);
  }

private:
  const Frees& frees_;

  std::ostream& os_;
  size_t& i_;

};

} // composite


template<typename VA>
void composite::Class<VA>::displayKey_v(std::ostream& os, size_t& i) const
{
  CALL_composite_worker( BOOST_PP_CAT(Display,KEY_name)<VA>(frees_,os,i) );
}


#undef KEY_name

#define NumberClassName BOOST_PP_CAT(N,NumberName)
#define numberFuncName  BOOST_PP_CAT(n,NumberName)

template<typename VA>
class composite::Class<VA>::NumberClassName
{
public:
  NumberClassName(const Frees& frees, size_t& num) : frees_(frees), num_(num) {}

  typedef size_t result_type;
  
  template<typename Act>
  size_t operator()(size_t s, const Act& act)
  {
    return s+act.numberFuncName();
  }

  template<typename T>
  void operator()(T) const
  {
    num_+=frees_(T::value).numberFuncName();
  }

private:
  const Frees& frees_;

  size_t& num_;

};



template<typename VA>
size_t composite::Class<VA>::BOOST_PP_CAT(numberFuncName,_v)() const
{
  size_t res=0;
  const NumberClassName helper(frees_,res);
  mpl::for_each<Ordinals>(helper);
  res=boost::fusion::fold(acts_,res,helper);
  return res;
}


#undef numberFuncName
#undef NumberClassName

#undef NumberName
#undef Class
#undef A_or_L

#undef COMPOSITE_AVERAGED_LIOUVILLEAN_ITERATING

#endif // COMPOSITE_AVERAGED_LIOUVILLEAN_ITERATING
