// -*- C++ -*-
#ifndef   _COMPOSITE_SYSTEM_IMPL_INCLUDED
#define   _COMPOSITE_SYSTEM_IMPL_INCLUDED

#include "Composite.h"

#include "Interaction.h"

#include "impl/LazyDensityOperatorSliceIterator.tcc"

#include "Exception.h"

#include "Algorithm.h"
#include "Range.h"

#include<boost/fusion/algorithm/iteration/for_each.hpp>
#include<boost/fusion/algorithm/iteration/fold.hpp>

#include<boost/mpl/for_each.hpp>

#include<boost/bind.hpp>

#include<boost/lambda/lambda.hpp>
namespace bll=boost::lambda;

#include<algorithm>
#include<list>


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


template<typename VA> template<typename H>
void Composite<VA>::worker(const H& helper) const
{
  mpl::for_each<Ordinals>(helper);
  boost::fusion::for_each(acts_,helper);

  // IMPORTANT NOTICE!!! Both mpl::for_each and
  // boost::fusion::for_each stores the helper BY VALUE. So, what
  // happens is that the initial helper is COPIED to both, so the
  // changes made during mpl::for_each in those members of the helper
  // which are stored BY VALUE in H, will be "undone" at startup of
  // boost::fusion::for_each. Ergo, there MUST NOT BE such
  // members. (There can be const members stored by value though, of
  // course.)
  //
  // How about static members?

}


///////////////
//
// Construction
//
///////////////


namespace composite {


using structure::SubSystemFree;


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

#define RETURN_type typename Composite<VA>::Frees

template<typename VA>
const RETURN_type
Composite<VA>::fillFrees(const VA& acts)
{
  RETURN_type res; res=structure::SubSystemFree();
  boost::fusion::for_each(acts,composite::FillFrees<RANK>(res));
  return res;
}

#undef  RETURN_type
#define RETURN_type typename Composite<VA>::Dimensions

template<typename VA>
const RETURN_type
Composite<VA>::fillDimensions(const Frees& frees)
{
  RETURN_type res;
  mpl::for_each<Ordinals>(composite::FillDimensions<RANK>(frees,res));
  return res;
}

#undef  RETURN_type



/////////////
//
// Parameters
//
/////////////



template<typename VA>
class Composite<VA>::DisplayParameters
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
void Composite<VA>::displayParameters(std::ostream& os) const
{
  os<<"# Composite\n# Dimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<std::endl;
  worker(DisplayParameters(frees_,os));
}

//////////////
//
// Frequencies
//
//////////////


namespace composite {

bool compFrees(const SubSystemFree& ssf1, const SubSystemFree& ssf2)
{
  return ssf1.get()->highestFrequency() < ssf2.get()->highestFrequency();
}

} // composite


template<typename VA>
double Composite<VA>::highestFrequency() const
{
  return boost::max_element(frees_,composite::compFrees)->get()->highestFrequency();
  // NEEDS_WORK add the interactions here
}


////////
//
// Exact
//
////////


template<typename VA>
class Composite<VA>::IsUnitary
{
public:
  IsUnitary(const Frees& frees, bool& isIt) : frees_(frees), isIt_(isIt) {}

  typedef bool result_type;
  
  template<typename Act>
  bool operator()(bool s, const Act& act)
  {
    return s && structure::ExactCommon::isUnitary(act.getEx());
  }

  template<typename T>
  void operator()(T) const
  {
    isIt_&=structure::ExactCommon::isUnitary(frees_(T::value).getEx());
  }

private:
  const Frees& frees_;

  bool& isIt_;

};



template<typename VA>
bool Composite<VA>::isUnitary() const
{
  bool res=true;
  IsUnitary helper(frees_,res);
  mpl::for_each<Ordinals>(helper);
  if (res) res=boost::fusion::fold(acts_,res,helper);
  return res;
}



template<typename VA>
class Composite<VA>::ActWithU
{
public:
  ActWithU(const Frees& frees, double dtdid, StateVectorLow& psi) : frees_(frees), dtdid_(dtdid), psi_(psi) {}

  template<typename Vec, typename Ex>
  void help(const Ex*const ex, Vec v) const
  {
    using namespace blitzplusplus::basi;

    if (ex) boost::for_each(fullRange(psi_,v),bind(&Ex::actWithU,ex,dtdid_,_1));
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getEx(),act);
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getEx(),tmptools::Vector<IDX>());
  }


private:
  const Frees& frees_;

  const double dtdid_;
  StateVectorLow& psi_; 

};


template<typename VA>
void Composite<VA>::actWithU(double dtdid, StateVectorLow& psi) const
{
  worker(ActWithU(frees_,dtdid,psi));
}


//////////////
//
// Hamiltonian
//
//////////////

template<typename VA>
class Composite<VA>::Hamiltonian
{
public:
  Hamiltonian(const Frees& frees, double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) 
    : frees_(frees), t_(t), psi_(psi), dpsidt_(dpsidt), tIntPic0_(tIntPic0) {}

  template<typename Vec, typename Ha>
  void help(const Ha*const ha, Vec v) const
  {
    using namespace blitzplusplus;

    if (ha) 
      cpputils::for_each(basi::fullRange(psi_,v),basi::begin(dpsidt_,v),
		      bind(&Ha::addContribution,ha,t_,_1,_2,tIntPic0_)); 
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getHa(),act);
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getHa(),tmptools::Vector<IDX>());
  }

private:
  const Frees& frees_;

  const double t_;
  const StateVectorLow& psi_; 
  StateVectorLow& dpsidt_; 
  const double tIntPic0_;

};


template<typename VA>
void Composite<VA>::addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  worker(Hamiltonian(frees_,t,psi,dpsidt,tIntPic0));
}


//////////////
//
// Liouvillean
//
//////////////


template<typename VA>
class Composite<VA>::NJumps
{
public:
  NJumps(const Frees& frees, size_t& num) : frees_(frees), num_(num) {}

  typedef size_t result_type;
  
  template<typename Act>
  size_t operator()(size_t s, const Act& act)
  {
    return s+structure::LiouvilleanCommon::nJumps(act.getLi());
  }

  template<typename T>
  void operator()(T) const
  {
    num_+=structure::LiouvilleanCommon::nJumps(frees_(T::value).getLi());
  }

private:
  const Frees& frees_;

  size_t& num_;

};


template<typename VA>
size_t Composite<VA>::nJumps() const
{
  size_t res=0;
  NJumps helper(frees_,res);
  mpl::for_each<Ordinals>(helper);
  res=boost::fusion::fold(acts_,res,helper);
  return res;
}




template<typename VA>
class Composite<VA>::Probas
{
public:

  typedef std::list<Probabilities> SeqProbabilities;
  typedef typename SeqProbabilities::iterator SAI;

  Probas(const Frees& frees, double t, const LazyDensityOperator& ldo, SAI& iter) 
    : frees_(frees), t_(t), ldo_(ldo), iter_(iter) {}

  template<typename Vec, typename Li>
  void help(const Li*const li, Vec v) const
  {
    iter_++->reference(quantumdata::partialTrace(ldo_,bind(&Li::probabilities,t_,_1,li,structure::theStaticOne),v,defaultArray));    
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getLi(),act);
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getLi(),tmptools::Vector<IDX>());
  }

private:
  const Frees& frees_;

  const double t_;
  const LazyDensityOperator& ldo_;
  
  SAI& iter_;

};


template<typename VA>
const typename Composite<VA>::Probabilities
Composite<VA>::probabilities(double t, const LazyDensityOperator& ldo) const
{
  using namespace std;
  list<Probabilities> seqProbabilities(RANK+mpl::size<VA>::value);

  {
    typename list<Probabilities>::iterator iter(seqProbabilities.begin());

    worker(Probas(frees_,t,ldo,iter));
  }

  Probabilities res(nJumps()); res=0;
  return cpputils::concatenate(seqProbabilities,res);

}



template<typename VA>
class Composite<VA>::ActWithJ
{
public:
  ActWithJ(const Frees& frees, double t, StateVectorLow& psi, size_t& ordoJump, bool& flag) : frees_(frees), t_(t), psi_(psi), ordoJump_(ordoJump), flag_(flag) {}

  template<typename Vec, typename Li>
  void help(const Li*const li, Vec v) const
  {
    using namespace blitzplusplus::basi;

    if (!flag_ && li) {
      size_t n=Li::nJumps(li);
      if (ordoJump_<n) {
	boost::for_each(fullRange(psi_,v),bind(&Li::actWithJ,li,t_,_1,ordoJump_));
	flag_=true;
      }
      ordoJump_-=n;  
    }
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getLi(),act);
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getLi(),tmptools::Vector<IDX>());
  }

private:
  const Frees& frees_;

  const double t_;
  StateVectorLow& psi_; 
  size_t& ordoJump_;
  bool& flag_;

};


template<typename VA>
void Composite<VA>::actWithJ(double t, StateVectorLow& psi, size_t ordoJump) const
{
  bool flag=false;
  worker(ActWithJ(frees_,t,psi,ordoJump,flag));
}



///////////
//
// Averaged
//
///////////


template<typename VA>
class Composite<VA>::DisplayKey
{
public:
  DisplayKey(const Frees& frees, std::ostream& os, size_t& i) 
    : frees_(frees), os_(os), i_(i) {}

  template<typename Act>
  void operator()(const Act& act) const
  {
    structure::AveragedCommon::displayKey(os_,i_,act.getAv());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    structure::AveragedCommon::displayKey(os_,i_,frees_(IDX).getAv());
  }

private:
  const Frees& frees_;

  std::ostream& os_;
  size_t& i_;

};


template<typename VA>
void Composite<VA>::displayKey(std::ostream& os, size_t& i) const
{
  worker(DisplayKey(frees_,os,i));
}



template<typename VA>
class Composite<VA>::NAvr
{
public:
  NAvr(const Frees& frees, size_t& num) : frees_(frees), num_(num) {}

  typedef size_t result_type;
  
  template<typename Act>
  size_t operator()(size_t s, const Act& act)
  {
    return s+structure::AveragedCommon::nAvr(act.getAv());
  }

  template<typename T>
  void operator()(T) const
  {
    num_+=structure::AveragedCommon::nAvr(frees_(T::value).getAv());
  }

private:
  const Frees& frees_;

  size_t& num_;

};



template<typename VA>
size_t Composite<VA>::nAvr() const
{
  size_t res=0;
  NAvr helper(frees_,res);
  mpl::for_each<Ordinals>(helper);
  res=boost::fusion::fold(acts_,res,helper);
  return res;
}




template<typename VA>
class Composite<VA>::Average
{
public:

  typedef std::list<Averages> SeqAverages;
  typedef typename SeqAverages::iterator SAI;

  Average(const Frees& frees, double t, const LazyDensityOperator& ldo, SAI& iter) 
    : frees_(frees), t_(t), ldo_(ldo), iter_(iter) {}

  template<typename Vec, typename Av>
  void help(const Av*const av, Vec v) const
  {
    iter_++->reference(quantumdata::partialTrace(ldo_,bind(&Av::average,t_,_1,av,structure::theStaticOne),v,defaultArray));
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getAv(),act);
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getAv(),tmptools::Vector<IDX>());
  }

private:
  const Frees& frees_;

  const double t_;
  const LazyDensityOperator& ldo_;
  
  SAI& iter_;

};


template<typename VA>
const typename Composite<VA>::Averages
Composite<VA>::average(double t, const LazyDensityOperator& ldo) const
{
  using namespace std;
  list<Averages> seqAverages(RANK+mpl::size<VA>::value);

  {
    typename list<Averages>::iterator iter(seqAverages.begin());

    worker(Average(frees_,t,ldo,iter));
  }

  Averages res(nAvr()); res=0;
  return cpputils::concatenate(seqAverages,res);

}




template<typename VA>
class Composite<VA>::Process
{
public:
  Process(const Frees& frees, const Averages& avr, ptrdiff_t& l, ptrdiff_t& u) 
    : frees_(frees), avr_(avr), l_(l), u_(u) {}


  template<typename Av>
  void help(const Av*const av) const
  {
    using blitz::Range;
    if ((u_=l_+Av::nAvr(av))>l_) {
      Averages temp(avr_(Range(l_+1,u_)));
      Av::process(temp,av);
    }
    std::swap(l_,u_);
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getAv());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getAv());
  }

private:
  const Frees& frees_;

  const Averages& avr_;

  ptrdiff_t &l_, &u_;

};


template<typename VA>
void
Composite<VA>::process(Averages& avr) const
{
  ptrdiff_t l=-1, u=0;
  worker(Process(frees_,avr,l,u));
}



template<typename VA>
class Composite<VA>::Display
{
public:
  Display(const Frees& frees, const Averages& avr, std::ostream& os, int precision, ptrdiff_t& l, ptrdiff_t& u) 
    : frees_(frees), avr_(avr), os_(os), precision_(precision), l_(l), u_(u) {}

  template<typename Av>
  void help(const Av*const av) const
  {
    using blitz::Range;
    if ((u_=l_+Av::nAvr(av))>l_) av->display(avr_(Range(l_+1,u_)),os_,precision_);
    std::swap(l_,u_);
  }

  template<typename Act>
  void operator()(const Act& act) const
  {
    help(act.getAv());
  }

  template<int IDX>
  void operator()(mpl::integral_c<int,IDX>) const
  {
    help(frees_(IDX).getAv());
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
Composite<VA>::display(const Averages& avr, std::ostream& os, int precision) const
{
  ptrdiff_t l=-1, u=0;
  worker(Display(frees_,avr,os,precision,l,u));
}



#endif // _COMPOSITE_SYSTEM_IMPL_INCLUDED
