// -*- C++ -*-
#ifndef   ELEMENTS_FREES_IMPL_MULTILEVEL_TCC_INCLUDED
#define   ELEMENTS_FREES_IMPL_MULTILEVEL_TCC_INCLUDED

#include "MultiLevel_.h"

#include<boost/assign/list_inserter.hpp>

#include<boost/lambda/bind.hpp>
#include<boost/lambda/lambda.hpp>
namespace bll=boost::lambda;

#include<boost/iterator/transform_iterator.hpp>

#include<boost/fusion/sequence/intrinsic/at_c.hpp>
// #include<boost/fusion/container/generation/make_list.hpp>
#include<boost/fusion/container/generation/make_vector.hpp>
#include<boost/fusion/algorithm/iteration/for_each.hpp>

#include<boost/fusion/mpl/at.hpp>

#include<boost/mpl/for_each.hpp>

#include<algorithm>
#include<numeric>


namespace multilevel {


// using boost::fusion::make_list;
using boost::fusion::make_vector;


namespace result_of {

// using boost::fusion::result_of::make_list;
using boost::fusion::result_of::make_vector;

} // result_of


////////
//
// Exact
//
////////


template<int NL>
void
Exact<NL>::updateU(double dtdid) const
{
  using namespace std;
  dcomp (*dcompExp)(const dcomp&)=&exp;
  boost::transform(zIs_,getFactors().begin(),bll::bind(dcompExp,-dtdid*bll::_1));
}


template<int NL>
bool 
Exact<NL>::isUnitary_v() const
{
  using namespace std  ;
  // using namespace boost::lambda;
  // using namespace lambda;

  return !accumulate(
		     boost::make_transform_iterator(zIs_.begin(),hasRealPart),
		     boost::make_transform_iterator(zIs_.end  (),hasRealPart),
		     false,
		     bll::_1 || bll::_2 //logical_or<bool>()
		     );
}


/*
template<int NL, int N1, int N2>
const quantumoperator::Frequencies<1>
shift(const MultiLevelBase<NL>& mlb, const Frequencies& freq)
{
  else return freq;

}
*/

//////////////
//
// Hamiltonian
//
//////////////


namespace details {


template<int NL>
struct ElementaryLevel
{
  typedef typename LevelsMF<NL>::type Levels;

  ElementaryLevel(const StateVectorLow& psi, StateVectorLow& dpsidt, const Levels& zs) : psi_(psi), dpsidt_(dpsidt), zs_(zs) {}

  template<typename T>
  void operator()(T) const
  {
    dpsidt_(T::value)+=-zs_(T::value)*psi_(T::value);
  }

private:
  const StateVectorLow& psi_;
  StateVectorLow& dpsidt_;
  const Levels& zs_;
  
};


template<int NL>
struct ElementaryPumping
{

  ElementaryPumping(const StateVectorLow& psi, StateVectorLow& dpsidt) : psi_(psi), dpsidt_(dpsidt) {}

  template<typename P>
  void operator()(const P& pump) const
  {
    typename P::template SanityCheck<0,NL-1>(); 
    // A temporary variable to instantiate the SanityCheck member template
    dpsidt_(P::first )+=conj(pump.get())*psi_(P::second);
    dpsidt_(P::second)-=     pump.get() *psi_(P::first );
  }

private:
  const StateVectorLow& psi_;
  StateVectorLow& dpsidt_;
  
};


} // details



template<int NL, typename VP>
void
HamiltonianSch<NL,VP>::addContribution_v(const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  using namespace details;
  mpl::for_each<tmptools::Ordinals<NL> >(ElementaryLevel<NL>(psi,dpsidt,zSchs_));
  for_each(etas_,ElementaryPumping<NL>(psi,dpsidt));
}


//////////////
//
// Liouvillean
//
//////////////


template<int NL, typename VL> template<int N>
void 
Liouvillean<NL,VL>::jumpStrategy(StateVectorLow& psi) const
{
  typedef typename mpl::at_c<VL, N>::type Decay;
  typename Decay::template SanityCheck<0,NL-1>();
  dcomp temp(sqrt(2.*boost::fusion::at_c<N>(gammas_).get())*psi(Decay::second));
  psi=0;
  psi(Decay::first)=temp;

}


template<int NL, typename VL> template<int N>
double
Liouvillean<NL,VL>::jumpProbabilityStrategy(const LazyDensityOperator& matrix) const
{
  typedef typename mpl::at_c<VL, N>::type Decay;
  return 2.*boost::fusion::at_c<N>(gammas_).get()*matrix(Decay::second);
}



template<int NL, typename VL>
class Liouvillean<NL,VL>::JS_helper
{
public:
  JS_helper(JumpStrategies& jumpStrategies, const Liouvillean* liouvillean)
    : jumpStrategies_(jumpStrategies), liouvillean_(liouvillean) {}

  template<typename T>
  void operator()(T) const
  {
    jumpStrategies_(T::value)=bind(&Liouvillean::template jumpStrategy<T::value>,
				   liouvillean_,
				   _1);
  }

private:
  JumpStrategies& jumpStrategies_;
  const Liouvillean* liouvillean_;
  
};


template<int NL, typename VL>
class Liouvillean<NL,VL>::JPS_helper
{
public:
  JPS_helper(JumpProbabilityStrategies& jumpProbabilityStrategies, const Liouvillean* liouvillean)
    : jumpProbabilityStrategies_(jumpProbabilityStrategies), liouvillean_(liouvillean) {}

  template<typename T>
  void operator()(T) const
  {
    jumpProbabilityStrategies_(T::value)=bind(&Liouvillean::template jumpProbabilityStrategy<T::value>,
					      liouvillean_,
					      _1);
  }

private:
  JumpProbabilityStrategies& jumpProbabilityStrategies_;
  const Liouvillean* liouvillean_;
  
};



template<int NL, typename VL>
const typename Liouvillean<NL,VL>::JumpStrategies
Liouvillean<NL,VL>::fillJS() const
{
  typename Liouvillean::JumpStrategies res;
  mpl::for_each<tmptools::Ordinals<NLT> >(JS_helper(res,this));
  return res;
}


template<int NL, typename VL>
const typename Liouvillean<NL,VL>::JumpProbabilityStrategies
Liouvillean<NL,VL>::fillJPS() const
{
  typename Liouvillean::JumpProbabilityStrategies res;
  mpl::for_each<tmptools::Ordinals<NLT> >(JPS_helper(res,this));
  return res;
}



template<int NL, typename VL>
class Liouvillean<NL,VL>::KeyHelper
{
public:
  KeyHelper(KeyLabels& keyLabels) : keyLabels_(keyLabels) {}

  template<typename T>
  void operator()(T)
  {
    using namespace std;
    stringstream slate(stringstream::out);
    slate<<"Jump "<<T::second<<" -> "<<T::first;
    keyLabels_.push_back(slate.str());
  }

private:
  KeyLabels& keyLabels_;

};


template<int NL, typename VL>
const typename Liouvillean<NL,VL>::KeyLabels
Liouvillean<NL,VL>::fillKeyLabels()
{
  typename Liouvillean::KeyLabels res;
  mpl::for_each<VL>(KeyHelper(res));
  return res;
}

///////////
//
// 
//
///////////



//////////
//
// Helpers
//
//////////

namespace details {


template<int NL>
struct ElementaryDL
{

  typedef blitz::TinyVector<dcomp,NL> Levels;

  ElementaryDL(Levels& levels) : levels_(levels) {}

  template<typename P>
  void operator()(const P& gamma) const
  {
    levels_(P::second)+=gamma.get();
  }

private:
  Levels& levels_;
  
};


} // details


#define RETURN_type blitz::TinyVector<dcomp,NL>

template<int NL, typename VL>
const RETURN_type
decayingLevels(const blitz::TinyVector<double,NL>& deltas, const VL& gammas)
{
  RETURN_type res(deltas);
  res*=-DCOMP_I;
  for_each(gammas,details::ElementaryDL<NL>(res));
  return res;
}

#undef  RETURN_type


template<int NL>
const structure::DynamicsBase::RealFreqs
filterReal(const blitz::TinyVector<dcomp,NL>& levels)
{
  using namespace std;
  structure::DynamicsBase::RealFreqs res;
  for (int i=0; i<NL; i++)
    if (!hasRealPart(levels(i))) {
      stringstream tag(stringstream::out);
      tag<<"delta"<<i;
      res.push_back(boost::make_tuple(tag.str(),-imag(levels(i)),1.));
    }
  return res;
}


struct ElementaryComplexFreqs
{
  typedef structure::DynamicsBase::ComplexFreqs ComplexFreqs;

  ElementaryComplexFreqs(ComplexFreqs& cf, const std::string& label) : cf_(cf), label_(label) {}

  template<typename P>
  void operator()(const P& eta) const
  {
    using namespace std;
    using boost::assign::insert;

    stringstream tag(stringstream::out);
    tag<<label_<<P::first<<P::second;
    cf_.push_back(boost::make_tuple(tag.str(),eta.get(),1.));
  }

private:
  ComplexFreqs& cf_;
  const std::string& label_;

};



#define RETURN_type structure::DynamicsBase::ComplexFreqs


template<int NL, typename VP>
const RETURN_type
complexFreqs(const blitz::TinyVector<dcomp,NL>& levels, const VP& etas)
{
  using namespace std;
  using boost::assign::insert;
  RETURN_type res;
  for (int i=0; i<NL; i++)
    if (hasRealPart(levels(i))) {
      stringstream tag(stringstream::out);
      tag<<"(gamma"<<i<<",-delta"<<i<<")";
      res.push_back(boost::make_tuple(tag.str(),levels(i),1.));
    }
  for_each(etas,ElementaryComplexFreqs(res,"eta"));
  return res;
}


#undef  RETURN_type


} // multilevel


////////////////
//
// Highest level
//
////////////////


template<int NL, typename VP, typename VL, typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedLossyMultiLevelSch<NL,VP,VL,AveragingType>::PumpedLossyMultiLevelSch(const RealLevels& deltas, const VP& etas, const VL& gammas, const AveragingConstructorParameters&... a)
  : Hamiltonian(multilevel::decayingLevels(deltas,gammas),etas),
    Liouvillean(gammas),
    Base(multilevel::filterReal(get_zSchs()),multilevel::complexFreqs(get_zSchs(),etas)),
    AveragingType(a...)
{
  getParsStream()<<"# Schroedinger picture.\n";
}


#endif // ELEMENTS_FREES_IMPL_MULTILEVEL_TCC_INCLUDED



/*
namespace multilevel {


template<int NL, typename VP>
inline
const HamiltonianIP<NL,VP>
makeHamiltonianIP(const blitz::TinyVector<dcomp,NL>& zSchs, const blitz::TinyVector<dcomp,NL> zIs, const VP& etas
		  )
{
  return HamiltonianIP<NL,VP>(zSchs,zIs,etas);
}



template<int NL, typename VP>
inline
const HamiltonianSch<NL,VP>
makeHamiltonianSch(const blitz::TinyVector<dcomp,NL>& zSchs, const VP& etas)
{
  return HamiltonianSch<NL,VP>(zSchs,etas);
}


template<int NL, typename VL>
inline
const Liouvillean<NL,VL>
makeLiouvillean(const VL& gammas)
{
  return Liouvillean<NL,VL>(gammas);
}


void foo()
{
  using namespace boost::fusion;

  makeHamiltonianIP(
		    LevelsMF<4>::type(dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2)),
		    LevelsMF<4>::type(dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2)),
		    make_list(
			      Pump<3,2>(dcomp(1.,-.2)),
			      Pump<3,1>(dcomp(1.,-.2))
			      )
		    );

  makeHamiltonianSch(
		     LevelsMF<4>::type(dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2)),
		     make_list(
			       Pump<3,2>(dcomp(1.,-.2)),
			       Pump<3,1>(dcomp(1.,-.2))
			       )
		     );

  makeLiouvillean(
		  make_list(
			    Decay<3,2>(-.2),
			    Decay<3,1>(1.)
			    )
		  );

}


} // multilevel


void foo()
{
  using namespace boost::fusion;
  using namespace multilevel;

  makePumpedLossyMultiLevelSch(
			       LevelsMF<4>::type(dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2),dcomp(1.,-.2)),
			       make_list(
					 Pump<3,2>(dcomp(1.,-.2)),
					 Pump<3,1>(dcomp(1.,-.2))
					 ),
			       make_list(
					 Decay<3,2>(-.2),
					 Decay<3,1>(1.)
					 ),
			       multilevel::EmptyBase<4>()
			       );

}

  
*/
