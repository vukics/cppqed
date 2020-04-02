/// \briefFile{Defines free elements of the \ref multilevelbundle "MultiLevel bundle" }
// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED
#define   CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED

#include "MultiLevel_Fwd.h"

#include "ParsMultiLevel.h"

#include "AveragingUtils.tcc"

#include "ElementAveraged.h"
#include "ElementLiouvillean.h"
#include "Free.h"
#include "FreeExact.h"
#include "Hamiltonian.h"

#include "primitive.hpp"

#include <boost/fusion/mpl/size.hpp>



/// Contains helpers for the \ref multilevelbundle "MultiLevel bundle"
namespace multilevel {


const std::string keyTitle="MultiLevel";


using namespace structure::freesystem; using structure::NoTime;



////////
//
// Exact
//
////////

struct MultiLevelExactNotImplementedException : public cpputils::Exception {};

/** \todo MultiLevel Exact not fully implemented, cf. the late 'shift' function. Maybe it's easier with something similar to Tridiagonal::freqs_ in Sigma */
template<int NL>
class Exact : public structure::FreeExact<false>
{
public:
  typedef ComplexPerLevel<NL> L;

  Exact(const L& zIs) : FreeExact<false>(NL), zIs_(zIs) {}

  const L& get_zIs() const {return zIs_;}

private:
  void updateU(double) const;

  bool applicableInMaster_v() const;

  const L zIs_;

};



//////////////
//
// Hamiltonian
//
//////////////


using cpputils::primitive;

/// Class representing an elementary pump term (an \f$\eta_{ij}\f$ \ref multilevelactualHamiltonian "here") with a compile-time pair \f$i,j\f$ and a runtime complex value
template<int I, int J>
class Pump : public primitive<dcomp>, public tmptools::pair_c<I,J>
{
public:
  using primitive<dcomp>::primitive;

};


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
  void addContribution_v(NoTime, const StateVectorLow&, StateVectorLow&) const;


  const L zSchs_;

  const VP etas_;

};


//////////////
//
// Liouvillean
//
//////////////


/// Class representing an elementary decay term (a \f$\gamma_{ij}\f$ \ref multilevelelements "here") with a compile-time pair \f$i,j\f$ and a runtime real value
template<int I, int J>
class Decay : public primitive<double>, public tmptools::pair_c<I,J>
{
public:
  using primitive<double>::primitive;

};


template<int NL, typename VL, bool IS_DIFFUSIVE, int ORDO>
class LiouvilleanRadiative : public LiouvilleanRadiative<NL,VL,IS_DIFFUSIVE,ORDO-1>
{
private:
  typedef LiouvilleanRadiative<NL,VL,IS_DIFFUSIVE,ORDO-1> Base;
  
protected:
  using Base::Base; using Base::NRT; using Base::gammas_;
  
  typedef typename Base::template LindbladNo<ORDO> LindbladOrdo; // shouldn’t be named LindbladNo, as then it would be confused with the LindbladNo in the direct base, 
                                                                 // which is not a template
  
private:
  void doActWithJ(NoTime, StateVectorLow&, LindbladOrdo) const override;
  double rate(NoTime, const LazyDensityOperator&, LindbladOrdo) const override;
  void doActWithSuperoperator(NoTime, const DensityOperatorLow&, DensityOperatorLow&, LindbladOrdo) const override;

};


namespace details {

template<int ORDO>
void doReallyActWithSuperoperator(const DensityOperatorLow&, DensityOperatorLow&, double) {} // trivial implementation for ORDO\neq0

} // details

template<int NL, typename VL, int ORDO>
class LiouvilleanDiffusive : public LiouvilleanDiffusive<NL,VL,ORDO-1>
{
private:
  typedef LiouvilleanDiffusive<NL,VL,ORDO-1> Base;
  
protected:
  using Base::Base; using Base::NRT; using Base::gamma_parallel_;
  
  typedef typename Base::template LindbladNo<NRT+ORDO> LindbladOrdo;
  
private:
  void doActWithJ(NoTime, StateVectorLow& psi, LindbladOrdo) const override {psi*=sqrt(2.*gamma_parallel_); psi(ORDO)*=-1.;}
  double rate(NoTime, const LazyDensityOperator&, LindbladOrdo) const override {return -1;}
  void doActWithSuperoperator(NoTime, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, LindbladOrdo) const override {details::doReallyActWithSuperoperator<ORDO>(rho,drhodt,gamma_parallel_);}

};


template<int NL, typename VL, bool IS_DIFFUSIVE>
class LiouvilleanBase;


// With ORDO=-1, these don’t correspond to valid jumps, they are here merely to store data:
#define BASE_class mpl::if_c<IS_DIFFUSIVE,LiouvilleanDiffusive<NL,VL,NL-1>,LiouvilleanBase<NL,VL,false> >::type
template<int NL, typename VL, bool IS_DIFFUSIVE>
class LiouvilleanRadiative<NL,VL,IS_DIFFUSIVE,-1> : public BASE_class
{
protected:
  LiouvilleanRadiative(const VL& gammas, double gamma_parallel) : BASE_class(gamma_parallel), gammas_(gammas) {}
#undef BASE_class

  const VL gammas_;
};


template<int NL, typename VL>
class LiouvilleanDiffusive<NL,VL,-1> : public LiouvilleanBase<NL,VL,true>
{
protected:
  LiouvilleanDiffusive(double gamma_parallel) : gamma_parallel_(gamma_parallel) {}
  
  const double gamma_parallel_;
  
};


#define BASE_class structure::ElementLiouvillean<1, mpl::size<VL>::value + tmptools::integral_if_c<IS_DIFFUSIVE,NL,0>::type::value >
template<int NL, typename VL, bool IS_DIFFUSIVE>
class LiouvilleanBase : public BASE_class // takes care of key composition
{
private:
  typedef BASE_class Base;
#undef BASE_class
  class KeyHelper;
  typedef typename Base::KeyLabels KeyLabels;
  static const KeyLabels fillKeyLabels();
  
protected:
  static const int NRT=mpl::size<VL>::value; // number of transitions with radiative loss

  LiouvilleanBase() : Base(keyTitle,fillKeyLabels()) {}
  LiouvilleanBase(double) : LiouvilleanBase() {}
  
};


#define BASE_class LiouvilleanRadiative<NL,VL,IS_DIFFUSIVE,mpl::size<VL>::value-1>
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


/// Implements a free multi-level system with driving and loss \see \ref multilevelbundle
/**
 * \tparam NL number of levels
 * \tparam VP fusion vector of multilevel::Pump types representing all the pump terms (the sequence of \f$\eta_{ij}\f$s \ref multilevelactualHamiltonian "here") in the system
 * \tparam VL fusion vector of multilevel::Decay types representing all the decays (the sequence of \f$\gamma_{ij}\f$s \ref multilevelelements "here") of the system
 *
 * \todo some static sanity checks of `VP` and `VL` in view of `NL` should be done
 */
template<int NL, typename VP, typename VL, bool IS_DIFFUSIVE, typename AveragingType>
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
  PumpedLossyMultiLevelSch(const RealPerLevel&, const VP&, const VL&, double gamma_parallel, AveragingConstructorParameters&&...);

};



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
