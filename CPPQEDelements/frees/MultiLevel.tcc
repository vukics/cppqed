// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_MULTILEVEL_TCC_INCLUDED
#define   CPPQEDELEMENTS_FREES_MULTILEVEL_TCC_INCLUDED

#include "MultiLevel_.h"

#include "ElementLiouvillean.tcc"

#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
namespace bll=boost::lambda;

#include <boost/range/adaptor/transformed.hpp>

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
// #include<boost/fusion/container/generation/make_list.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>

#include <boost/fusion/mpl/at.hpp>

#include <boost/mpl/for_each.hpp>


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
  dcomp (*dcompExp)(const dcomp&)=&exp;
  boost::transform(zIs_,getDiagonal().begin(),bll::bind(dcompExp,-dtdid*bll::_1));
}


template<int NL>
bool 
Exact<NL>::applicableInMaster_v() const
{
  return !accumulate(zIs_.begin() | boost::adaptors::transformed(hasRealPart),false,std::logical_or<bool>());
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
  typedef ComplexPerLevel<NL> L;

  ElementaryLevel(const StateVectorLow& psi, StateVectorLow& dpsidt, const L& zs) : psi_(psi), dpsidt_(dpsidt), zs_(zs) {}

  template<typename T>
  void operator()(T) const
  {
    dpsidt_(T::value)+=-zs_(T::value)*psi_(T::value);
  }

private:
  const StateVectorLow& psi_;
  StateVectorLow& dpsidt_;
  const L& zs_;
  
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
HamiltonianSch<NL,VP>::addContribution_v(NoTime, const StateVectorLow& psi, StateVectorLow& dpsidt) const
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


template<int NL, typename VL, bool IS_DIFFUSIVE, int ORDO>
void
LiouvilleanRadiative<NL,VL,IS_DIFFUSIVE,ORDO>::doActWithJ(NoTime, StateVectorLow& psi, LindbladOrdo) const
{
  typedef typename mpl::at_c<VL,ORDO>::type Decay;
  typename Decay::template SanityCheck<0,NL-1>();
  dcomp temp(sqrt(2.*boost::fusion::at_c<ORDO>(gammas_).get())*psi(Decay::second));
  psi=0;
  psi(Decay::first)=temp;
}


template<int NL, typename VL, bool IS_DIFFUSIVE, int ORDO>
double
LiouvilleanRadiative<NL,VL,IS_DIFFUSIVE,ORDO>::rate(NoTime, const LazyDensityOperator& matrix, LindbladOrdo) const
{
  typedef typename mpl::at_c<VL,ORDO>::type Decay;
  return 2.*boost::fusion::at_c<ORDO>(gammas_).get()*matrix(Decay::second);
}


template<int NL, typename VL, bool IS_DIFFUSIVE>
class LiouvilleanBase<NL,VL,IS_DIFFUSIVE>::KeyHelper
{
public:
  KeyHelper(KeyLabels& keyLabels) : keyLabels_(keyLabels) {}

  template<typename T>
  void operator()(T)
  {
    using namespace std;
    ostringstream slate;
    slate<<"Jump "<<T::second<<" -> "<<T::first;
    keyLabels_.push_back(slate.str());
  }

private:
  KeyLabels& keyLabels_;

};


template<int NL, typename VL, bool IS_DIFFUSIVE>
auto
LiouvilleanBase<NL,VL,IS_DIFFUSIVE>::fillKeyLabels() -> const KeyLabels
{
  KeyLabels res;
  mpl::for_each<VL>(KeyHelper(res));
  if (IS_DIFFUSIVE) for (int i=0; i<NL; ++i) { // For the moment simply implemented as a runtime loop (dirty but simple)
    ostringstream slate;
    slate<<"Phase flip for level "<<i;
    res.push_back(slate.str());
  }
  return res;
}

//////////
//
// Helpers
//
//////////

namespace details {


template<int NL>
struct ElementaryDL
{

  typedef ComplexPerLevel<NL> L;

  ElementaryDL(L& levels) : levels_(levels) {}

  template<typename P>
  void operator()(const P& gamma) const
  {
    levels_(P::second)+=gamma.get();
  }

private:
  L& levels_;
  
};


} // details


template<int NL, typename VL>
const auto
decayingLevels(const RealPerLevel<NL>& deltas, const VL& gammas)
{
  blitz::TinyVector<dcomp,NL> res(deltas);
  res*=-DCOMP_I;
  for_each(gammas,details::ElementaryDL<NL>(res));
  return res;
}


template<int NL>
const auto
filterReal(const ComplexPerLevel<NL>& levels, double gamma_parallel)
{
  using namespace std;
  structure::DynamicsBase::RealFreqs res;
  for (int i=0; i<NL; i++)
    if (!hasRealPart(levels(i))) {
      ostringstream tag;
      tag<<"delta"<<i;
      res.push_back(make_tuple(tag.str(),-imag(levels(i)),1.));
    }
  res.push_back(make_tuple("gamma_parallel",gamma_parallel,1.));
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

    ostringstream tag;
    tag<<label_<<P::first<<P::second;
    cf_.push_back(make_tuple(tag.str(),eta.get(),1.));
  }

private:
  ComplexFreqs& cf_;
  const std::string& label_;

};



template<int NL, typename VP>
const auto
complexFreqs(const ComplexPerLevel<NL>& levels, const VP& etas)
{
  using namespace std;

  structure::DynamicsBase::ComplexFreqs res;
  for (int i=0; i<NL; i++)
    if (hasRealPart(levels(i))) {
      ostringstream tag;
      tag<<"(gamma"<<i<<",-delta"<<i<<")";
      res.push_back(make_tuple(tag.str(),levels(i),1.));
    }
  for_each(etas,ElementaryComplexFreqs(res,"eta"));
  return res;
}



} // multilevel


////////////////
//
// Highest level
//
////////////////


template<int NL, typename VP, typename VL, bool IS_DIFFUSIVE, typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedLossyMultiLevelSch<NL,VP,VL,IS_DIFFUSIVE,AveragingType>::PumpedLossyMultiLevelSch(const RealPerLevel& deltas, const VP& etas, const VL& gammas, double gamma_parallel, AveragingConstructorParameters&&... a)
  : Hamiltonian(multilevel::decayingLevels(deltas,gammas),etas),
    Liouvillean(gammas,gamma_parallel),
    Base(multilevel::filterReal(this->get_zSchs(),gamma_parallel),multilevel::complexFreqs(this->get_zSchs(),etas)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  this->getParsStream()<<"Schroedinger picture.\n";
}


#endif // CPPQEDELEMENTS_FREES_MULTILEVEL_TCC_INCLUDED

