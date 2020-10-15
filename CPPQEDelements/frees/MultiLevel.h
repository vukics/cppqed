/// \briefFile{Defines free elements of the \ref multilevelbundle "MultiLevel bundle" }
// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED
#define   CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED

#include "ParsMultiLevel.h"

#include "AveragingUtils.tcc"

#include "ElementAveraged.h"
#include "ElementLiouvillean.h"
#include "Free.h"
#include "FreeExact.h"
#include "Hamiltonian.h"

#include "primitive.hpp"

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
// #include<boost/fusion/container/generation/make_list.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>

#include <boost/fusion/mpl/at.hpp>
#include <boost/fusion/mpl/size.hpp>

#include <boost/mpl/for_each.hpp>


/*

Convention is the following: an (l,m) pair in VP (or VL) below means a sigma_lm=|l><m|

This yields an 

H_elementaryPumping=i*(conj(eta)*sigma-eta*sigmadagger)

term in the Hamiltonian

If a pair (l,m) is in VL, this means a jump with jump operator
sigma_lm. This means that the rate of the jump is proportional
to 

<sigmadagger_lm*sigma_lm>=<|m><l|*|l><m|>=<|m><m|>=rho_mm

*/

/// Contains helpers for the \ref multilevelbundle "MultiLevel bundle"
namespace multilevel {


const std::string keyTitle="MultiLevel";


using namespace structure::freesystem; using structure::NoTime;

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

/** \todo MultiLevel Exact not fully implemented, cf. the late 'shift' function. Maybe it's easier with something similar to Tridiagonal::freqs_ in Sigma */
template<int NL>
class Exact : public structure::FreeExact<false>
{
public:
  typedef ComplexPerLevel<NL> L;

  Exact(const L& zIs) : FreeExact<false>(NL), zIs_(zIs) {}

  const L& get_zIs() const {return zIs_;}

private:
  void updateU(double dtdid) const {
    boost::transform(zIs_,getDiagonal().begin(),[=](const dcomp& zI) {return exp(-zI*dtdid);});
  }

  bool applicableInMaster_v() const {
    return !boost::accumulate(zIs_.begin(),false,[](bool init, const dcomp& zI) {return init || hasRealPart(zI);});
  }

  const L zIs_;

};



//////////////
//
// Hamiltonian
//
//////////////


using cpputils::primitive;


template<typename T, int I, int J>
class DynamicsPair : public primitive<T>, public tmptools::pair_c<I,J>
{
public:
  using primitive<T>::primitive;

};


/// Class representing an elementary pump term (an \f$\eta_{ij}\f$ \ref multilevelactualHamiltonian "here") with a compile-time pair \f$i,j\f$ and a runtime complex value
template<int I, int J> using Pump = DynamicsPair<dcomp,I,J>;


template<int NL, typename VP>
class HamiltonianIP 
  : public structure::Hamiltonian<1>,
    public Exact<NL>
{
public:
  static const int NPT=mpl::size<VP>::value; // number of pumped transitions

  typedef typename Exact<NL>::L L;

  HamiltonianIP(const L& zSchs, const L& zIs, const VP& etas)
    : Exact<NL>(zIs), zSchs_(zSchs), etas_(etas) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const;


  const L zSchs_;

  const VP etas_;

};



template<int NL, typename VP>
class HamiltonianSch 
  : public structure::HamiltonianTimeDependenceDispatched<1,structure::NO_TIME>
{
public:
  static const int NPT=mpl::size<VP>::value; // number of pumped transitions

  typedef ComplexPerLevel<NL> L;

  HamiltonianSch(const L& zSchs, const VP& etas) : zSchs_(zSchs), etas_(etas) {}

  const L& get_zSchs() const {return zSchs_;}

private:
  void addContribution_v(NoTime, const StateVectorLow& psi, StateVectorLow& dpsidt) const {
    mpl::for_each<tmptools::Ordinals<NL>>([&](auto arg) {const int ind=decltype(arg)::value; dpsidt(ind)+=-zSchs_(ind)*psi(ind);});
    for_each(etas_,[&](auto pump) {
      using P=decltype(pump);
      typename P::template SanityCheck<0,NL-1>(); // A temporary variable to instantiate the SanityCheck member template
      dpsidt(P::first )+=conj(pump.get())*psi(P::second);
      dpsidt(P::second)-=     pump.get() *psi(P::first );
    });
  }

  const L zSchs_;

  const VP etas_;

};


//////////////
//
// Liouvillean
//
//////////////


/// Class representing an elementary decay term (a \f$\gamma_{ij}\f$ \ref multilevelelements "here") with a compile-time pair \f$i,j\f$ and a runtime real value
template<int I, int J> using Decay = DynamicsPair<double,I,J>;


template<int NL, typename VL, int ORDO>
class LiouvilleanRadiative;


// With ORDO=-1, this doesn’t correspond to valid jumps, its function is to store data, fill keyLabels, etc.:
template<int NL, typename VL>
class LiouvilleanRadiative<NL,VL,-1> : public structure::ElementLiouvillean<1, mpl::size<VL>::value>
{
private:
  typedef structure::ElementLiouvillean<1, mpl::size<VL>::value> Base;

  typedef typename Base::KeyLabels KeyLabels;
  
protected:
  static const int NRT=mpl::size<VL>::value; // number of transitions with radiative loss

  LiouvilleanRadiative(const VL& gammas,
                       double = 0. ///< dummy to conform with LiouvilleanDiffusive
                      )
    : Base(keyTitle,[] {
        KeyLabels res;
        mpl::for_each<VL>([&](auto arg) {res.push_back("Jump "+std::to_string(decltype(arg)::second)+" -> "+std::to_string(decltype(arg)::first));});
        return res;
      } () /* this double parenthesis syntax is required to actually run the lambda */ ),
      gammas_(gammas) {}
  
  const VL gammas_;

};



template<int NL, typename VL, int ORDO>
class LiouvilleanRadiative : public LiouvilleanRadiative<NL,VL,ORDO-1>
{
private:
  typedef LiouvilleanRadiative<NL,VL,ORDO-1> Base;
  
protected:
  using Base::Base; using Base::NRT; using Base::gammas_;
  
private:
  void doActWithJ(NoTime, StateVectorLow& psi, typename Base::template LindbladNo<ORDO>) const override
  {
    typedef typename mpl::at_c<VL,ORDO>::type Decay;
    typename Decay::template SanityCheck<0,NL-1>();
    dcomp temp(sqrt(2.*boost::fusion::at_c<ORDO>(gammas_).get())*psi(Decay::second));
    psi=0;
    psi(Decay::first)=temp;
  }
  
  double rate(NoTime, const LazyDensityOperator& matrix, typename Base::template LindbladNo<ORDO>) const override
  {
    typedef typename mpl::at_c<VL,ORDO>::type Decay;
    return 2.*boost::fusion::at_c<ORDO>(gammas_).get()*matrix(Decay::second);    
  }
  
  void doActWithSuperoperator(NoTime, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, typename Base::template LindbladNo<ORDO>) const override
  {
    typedef typename mpl::at_c<VL,ORDO>::type Decay;
    drhodt(Decay::first,Decay::first)+=2.*boost::fusion::at_c<ORDO>(gammas_).get()*rho(Decay::second,Decay::second);
  }

};


// Just a convenience layer above ElementLiouvilleanDiffusive
#define BASE_class structure::ElementLiouvilleanDiffusive<1,LiouvilleanRadiative<NL,VL,mpl::size<VL>::value-1>>
template<int NL, typename VL>
class LiouvilleanDiffusive : public BASE_class
{
protected:
  LiouvilleanDiffusive(const VL& gammas, double gamma_parallel) : BASE_class(NL,gamma_parallel,gammas) {}
#undef BASE_class
  
};


#define BASE_class std::conditional_t<IS_DIFFUSIVE,LiouvilleanDiffusive<NL,VL>,LiouvilleanRadiative<NL,VL,mpl::size<VL>::value-1>>
template<int NL, typename VL, bool IS_DIFFUSIVE>
class Liouvillean : public BASE_class
{
private:
  typedef BASE_class Base;
#undef BASE_class
  
protected:
  using Base::Base;
};


} // multilevel


////////////////
//
// Highest level
//
////////////////


template<int NL>
class MultiLevelBase 
  : public structure::Free
{
public:
  typedef std::shared_ptr<const MultiLevelBase> Ptr;

  using structure::Free::getParsStream;

  MultiLevelBase(const RealFreqs& realFreqs=RealFreqs(), const ComplexFreqs& complexFreqs=ComplexFreqs())
    : structure::Free(NL,realFreqs,complexFreqs)
  {
    getParsStream()<<multilevel::keyTitle<<std::endl;
  }

  virtual ~MultiLevelBase() {}

};


namespace multilevel {

auto elementaryComplexFreqs(structure::DynamicsBase::ComplexFreqs& cf, const std::string& label)
{
  return [&](auto eta) {cf.push_back({label+std::to_string(decltype(eta)::first)+std::to_string(decltype(eta)::second),eta.get(),1.});} ;
}

} // multilevel



/// Implements a free multi-level system with driving and loss \see \ref multilevelbundle
/**
 * \tparam NL number of levels
 * \tparam VP fusion vector of multilevel::Pump types representing all the pump terms (the sequence of \f$\eta_{ij}\f$s \ref multilevelactualHamiltonian "here") in the system
 * \tparam VL fusion vector of multilevel::Decay types representing all the decays (the sequence of \f$\gamma_{ij}\f$s \ref multilevelelements "here") of the system
 *
 * \todo some static sanity checks of `VP` and `VL` in view of `NL` should be done
 */
template<int NL, typename VP, typename VL, bool IS_DIFFUSIVE, typename AveragingType=ReducedDensityOperator<1>>
class PumpedLossyMultiLevelSch 
  : public multilevel::HamiltonianSch<NL,VP>,
    public multilevel::Liouvillean<NL,VL,IS_DIFFUSIVE>,
    public MultiLevelBase<NL>,
    public AveragingType
  // The ordering becomes important here
{
public:
  typedef multilevel::ComplexPerLevel<NL> ComplexPerLevel;
  typedef multilevel::   RealPerLevel<NL>    RealPerLevel;

  typedef multilevel::HamiltonianSch<NL,VP>              Hamiltonian;
  typedef multilevel::Liouvillean   <NL,VL,IS_DIFFUSIVE> Liouvillean;
  typedef MultiLevelBase<NL> Base;

  typedef typename Base::Ptr Ptr;

  template<typename... AveragingConstructorParameters>
  PumpedLossyMultiLevelSch(const RealPerLevel& deltas, const VP& etas, const VL& gammas, double gamma_parallel, AveragingConstructorParameters&&... a)
    : Hamiltonian([&] {
        ComplexPerLevel res(deltas); res*=-DCOMP_I;
        for_each(gammas,[&](auto gamma) {res(decltype(gamma)::second)+=gamma.get();});
        return res;
      } (), etas),
      Liouvillean(gammas,gamma_parallel),
      Base([this,gamma_parallel] {
        structure::DynamicsBase::RealFreqs res;
        for (int i=0; i<NL; i++) if (!hasRealPart(this->get_zSchs()(i))) res.push_back({"delta"+std::to_string(i),-imag(this->get_zSchs()(i)),1.});
        res.push_back({"gamma_parallel",gamma_parallel,1.});
        return res;
      } () , [this,&etas] {
        structure::DynamicsBase::ComplexFreqs res;
        for (int i=0; i<NL; i++) if (hasRealPart(this->get_zSchs()(i))) res.push_back({"(gamma"+std::to_string(i)+",-delta"+std::to_string(i)+")",this->get_zSchs()(i),1.});
        for_each(etas,multilevel::elementaryComplexFreqs(res,"eta"));
        return res;        
      } () ),
      AveragingType(std::forward<AveragingConstructorParameters>(a)...)
  {
    this->getParsStream()<<"Schroedinger picture.\n";
  }

};
// VP is a compile-time container of pairs, specifying which transitions are pumped. It should model a Boost.Fusion sequence,
// which stores the pairs for compile-time use and the pump Rabi frequencies for run-time use.



namespace multilevel {

#define RETURN_type typename MultiLevelBase<NL>::Ptr

/// Maker function for PumpedLossyMultiLevelSch
template<typename AveragingType, int NL, typename VP, typename VL, typename... AveragingConstructorParameters>
inline
RETURN_type
makePumpedLossySch(const RealPerLevel<NL>& deltas,
                   const VP& etas, const VL& gammas, double gamma_parallel,
                   AveragingConstructorParameters&&... a)
{
  if (gamma_parallel) return std::make_shared<PumpedLossyMultiLevelSch<NL,VP,VL,true ,AveragingType> >(deltas,etas,gammas,gamma_parallel,a...);
  else                return std::make_shared<PumpedLossyMultiLevelSch<NL,VP,VL,false,AveragingType> >(deltas,etas,gammas,gamma_parallel,a...);
}

/// \overload
template<typename AveragingType, int NL, typename VP, typename VL, typename... AveragingConstructorParameters>
inline
RETURN_type
makePumpedLossySch(const multilevel::ParsPumpedLossy<NL,VP,VL>& p, AveragingConstructorParameters&&... a)
{
  return makePumpedLossySch<AveragingType>(p.deltas,p.etas,p.gammas,p.gamma_parallel,a...);
}

/// \overload
template<int NL, typename VP, typename VL>
inline
RETURN_type
makePumpedLossySch(const RealPerLevel<NL>& deltas, const VP& etas, const VL& gammas, double gamma_parallel, const std::string& keyTitle="PumpedLossyMultiLevelSch", bool offDiagonals=false)
{
  return makePumpedLossySch<ReducedDensityOperator<1> >(deltas,etas,gammas,gamma_parallel,keyTitle,NL,offDiagonals);
}

/// \overload
template<int NL, typename VP, typename VL>
inline
RETURN_type
makePumpedLossySch(const multilevel::ParsPumpedLossy<NL,VP,VL>& p, const std::string& keyTitle="PumpedLossyMultiLevelSch", bool offDiagonals=false)
{
  return makePumpedLossySch<ReducedDensityOperator<1> >(p,keyTitle,NL,offDiagonals);
}


#undef  RETURN_type

} // multilevel


#endif // CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED
