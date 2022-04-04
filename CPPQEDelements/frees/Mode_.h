// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the \ref genericelementsfreesmode bundle (tackling the dynamics of a single harmonic-oscillator mode)}
#ifndef CPPQEDELEMENTS_FREES_MODE__H_INCLUDED
#define CPPQEDELEMENTS_FREES_MODE__H_INCLUDED

#include "ParsMode.h" // for the user interface

#include "QM_Picture.h"
#include "StateVector.h"

#include "ElementLiouvillean.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"


class ModeBase;

/// Contains helpers for the \ref genericelementsfreesmode bundle
namespace mode {


struct DoNotAverage {};


const std::string keyTitle="Mode";


using namespace structure::freesystem; using structure::NoTime;

typedef std::shared_ptr<const ModeBase> Ptr;

const Tridiagonal aop(size_t dim);
const Tridiagonal aop(Ptr);
const Tridiagonal nop(Ptr);

inline const Tridiagonal xop(Ptr mode) {return tridiagPlusHC(aop(mode))/sqrt(2.);}
// inline const Tridiagonal yop(mode::Ptr) {return ...}

/// Coherent state
/**
 * The implementation relies on mathutils::coherentElement, which works also for high Fock-state elements
 *
 * \note The user has to take care that `alpha` is not too large for the given `cutoff` (rule of thumb: `cutoff>|alpha|^2`)
 */
StateVector coherent(dcomp alpha, ///< amplitude
                     size_t cutoff ///< cutoff
                     );

StateVector fock(size_t n, size_t dim, double phase=0);

/// Dispatcher for initial condition
StateVector init(const Pars&);


template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const Pars           &, QM_Picture, AveragingConstructorParameters&&... );

template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsLossy      &, QM_Picture, AveragingConstructorParameters&&... );

template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsPumped     &, QM_Picture, AveragingConstructorParameters&&... );

template<typename AveragingType, typename... AveragingConstructorParameters>
const Ptr make(const ParsPumpedLossy&, QM_Picture, AveragingConstructorParameters&&... );


const Ptr make(const Pars           &, QM_Picture);
const Ptr make(const ParsLossy      &, QM_Picture);
const Ptr make(const ParsPumped     &, QM_Picture);
const Ptr make(const ParsPumpedLossy&, QM_Picture);


double photonNumber(const StateVectorLow&); 
// This can be implemented in an extremely elegant way using tensor notation from blitz++, but it is not enough, the one below is needed.
double photonNumber(const LazyDensityOperator&);


template<bool B>
inline std::ostream& isFiniteTempStream(std::ostream& os, double nTh)
{
  if constexpr (B) return os<<" Finite temperature.\nnTh="<<nTh<<std::endl;
  else return os<<std::endl;
}

template<bool B>
inline double finiteTemperatureHamiltonianDecay(const ParsLossy& p)
{
  if constexpr (B) return p.kappa*(2.*p.nTh+1);
  else return p.kappa;
}

const Tridiagonal::Diagonal mainDiagonal(dcomp z, double omegaKerr, size_t dim);

const Tridiagonal pumping(dcomp eta, size_t dim);



class Exact : public structure::FreeExact<false>
{
public:
  Exact(dcomp zI, double omegaKerr, size_t);

  dcomp get_zI() const {return zI_;}
  
  double get_omegaKerr() const {return omegaKerr_;}
  
private:
  void updateU(Time) const; ///< `Time` is structure::OneTime in this case

  bool applicableInMaster_v() const {return !bool(real(zI_));}

  const dcomp zI_;
  
  const double omegaKerr_;

};


namespace details { struct EmptyBase {}; }

template<bool IS_TIME_DEPENDENT>
class Hamiltonian 
  : public quantumoperator::TridiagonalHamiltonian<1,IS_TIME_DEPENDENT>,
    public std::conditional_t<IS_TIME_DEPENDENT,Exact,details::EmptyBase>
{
public:
  typedef quantumoperator::TridiagonalHamiltonian<1,IS_TIME_DEPENDENT> Base;

  template<bool B=IS_TIME_DEPENDENT> Hamiltonian(std::enable_if_t<B==true,dcomp> zSch, dcomp zI, dcomp eta, double omegaKerr, double omegaKerrAlter, size_t);
  template<bool B=IS_TIME_DEPENDENT> Hamiltonian(std::enable_if_t<B==false,dcomp> zSch, dcomp eta, double omegaKerr, double omegaKerrAlter, size_t);

protected:
  dcomp get_zSch() const {return zSch_;}
  dcomp get_zI  () const {return Exact::get_zI();}

  dcomp get_eta() const {return eta_;}
  
  double get_omegaKerr() const {return Exact::get_omegaKerr();}

private:
  const dcomp
    zSch_,
  // z_=kappa_-I*delta_ --- zSch_ is the part appearing in the tridiagonal Hamiltonian, zI_ is appearing in the frequencies.
    eta_;

  const size_t dim_;

};



template<bool TEMPERATURE, bool IS_ALTERNATIVE=false> class Liouvillean;


namespace details {

void aJump   (StateVectorLow&, double);
void aDagJump(StateVectorLow&, double);

void aSuperoperator   (const DensityOperatorLow&, DensityOperatorLow&, double);
void aDagSuperoperator(const DensityOperatorLow&, DensityOperatorLow&, double);

} // details


template<bool IS_ALTERNATIVE>
class Liouvillean<false,IS_ALTERNATIVE> 
  : public structure::ElementLiouvillean<1,1>
{
protected:
  Liouvillean(double kappa, double=0, const std::string& kT=keyTitle) : structure::ElementLiouvillean<1,1>(kT,"excitation loss"), kappa_(kappa) {}
  // the second dummy argument is there only to have the same form for the ctor as in the TEMPERATURE=true case

private:
  void   doActWithJ (NoTime, StateVectorLow& psi       ) const override {details::aJump(psi,kappa_);}
  
  double rate       (NoTime, const LazyDensityOperator& m) const override
  {
    if constexpr (IS_ALTERNATIVE) return -1.;
    else return 2*kappa_*photonNumber(m);
  }

  void doActWithSuperoperator(NoTime, const DensityOperatorLow& rho, DensityOperatorLow& drhodt) const override {details::aSuperoperator(rho,drhodt,kappa_);}
  
  const double kappa_;

};


namespace details {

class LiouvilleanFiniteTemperatureBase
{
protected:
  LiouvilleanFiniteTemperatureBase(double kappa, double nTh) : kappa_(kappa), nTh_(nTh) {}
  
  template<bool IA> double rate0(const LazyDensityOperator& m) const
  {
    if constexpr (IA) return -1.;
    else return 2.*kappa_*(nTh_+1)*photonNumber(m);
  }
  
  template<bool IA> double rate1(const LazyDensityOperator& m) const
  {
    if constexpr (IA) return -1.;
    else return 2.*kappa_*nTh_*(photonNumber(m)+m.trace());

  }
  
  const double kappa_, nTh_;

};

} // details


template<bool IS_ALTERNATIVE> 
class Liouvillean<true ,IS_ALTERNATIVE>
  : private details::LiouvilleanFiniteTemperatureBase,
    public structure::ElementLiouvillean<1,2>
{
protected:
  typedef structure::ElementLiouvillean<1,2> Base;

  Liouvillean(double kappa, double nTh, const std::string& kT=keyTitle)
    : details::LiouvilleanFiniteTemperatureBase(kappa,nTh), Base(kT,{"excitation loss","excitation absorption"}) {}
  
private:
  void doActWithJ(NoTime, StateVectorLow& psi, LindbladNo<0>) const override {details::   aJump(psi,kappa_*(nTh_+1));}
  void doActWithJ(NoTime, StateVectorLow& psi, LindbladNo<1>) const override {details::aDagJump(psi,kappa_* nTh_   );}
  
  double rate(NoTime, const LazyDensityOperator& matrix, LindbladNo<0>) const override {return rate0<IS_ALTERNATIVE>(matrix);}
  double rate(NoTime, const LazyDensityOperator& matrix, LindbladNo<1>) const override {return rate1<IS_ALTERNATIVE>(matrix);}

  void doActWithSuperoperator(NoTime, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, LindbladNo<0>) const override {details::   aSuperoperator(rho,drhodt,kappa_*(nTh_+1));}
  void doActWithSuperoperator(NoTime, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, LindbladNo<1>) const override {details::aDagSuperoperator(rho,drhodt,kappa_* nTh_   );}

};



// A basic, extensible Averaging class
class Averaged : public structure::ClonableElementAveraged<1>
{
public:
  typedef structure::ClonableElementAveraged<1> Base;

  Averaged(const KeyLabels& follow=KeyLabels(), const KeyLabels& precede=KeyLabels());

protected:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;
  void           process_v(        Averages&)                  const;

private:
  const ClonedPtr do_clone() const {return new Averaged(*this);}

};


class AveragedQuadratures : public Averaged
{
public:
  AveragedQuadratures(const KeyLabels& follow=KeyLabels(), const KeyLabels& precede=KeyLabels());

protected:
  const Averages average_v(NoTime, const LazyDensityOperator&) const;
  void           process_v(        Averages&)                  const;

private:
  const ClonedPtr do_clone() const {return new AveragedQuadratures(*this);}

};


template<typename Base=Averaged>
class AveragedMonitorCutoff : public Base
{
public:
  AveragedMonitorCutoff();

private:
  const Averages average_v(NoTime t, const LazyDensityOperator& matrix) const
  {
    auto averages{Base::average_v(t,matrix)}; // This is already of the correct size, since nAvr knows about the size updated by the derived class
    averages(averages.size()-1)=matrix(matrix.getDimension()-1);
    return averages;
  }
  
  void process_v(Averages&) const;

};


} // mode



////////////////
//
// Highest level
//
////////////////


class ModeBase 
  : public structure::Free
{
public:
  ModeBase(size_t dim,
           const RealFreqs& realFreqs=emptyRF, const ComplexFreqs& complexFreqs=emptyCF,
           const std::string& keyTitle=mode::keyTitle);

  ModeBase(size_t dim, const ComplexFreqs& complexFreqs, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,emptyRF,complexFreqs,keyTitle) {}
  ModeBase(size_t dim, RealFreqsInitializer rf, ComplexFreqsInitializer cf={}, const std::string& keyTitle=mode::keyTitle)
    : ModeBase(dim,RealFreqs(rf),ComplexFreqs(cf),keyTitle) {}
  ModeBase(size_t dim, ComplexFreqsInitializer cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,{},cf,keyTitle) {}
  ModeBase(size_t dim, RF rf, CF cf=CF(), const std::string& keyTitle=mode::keyTitle)
    : ModeBase(dim,RealFreqsInitializer{rf}, cf==CF() ? ComplexFreqsInitializer{} : ComplexFreqsInitializer{cf},keyTitle) {}
  ModeBase(size_t dim, CF cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,ComplexFreqsInitializer{cf},keyTitle) {}
  ModeBase(size_t dim, RealFreqsInitializer rf, CF cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,rf,{cf},keyTitle) {}
  ModeBase(size_t dim, RF rf, ComplexFreqsInitializer cf, const std::string& keyTitle=mode::keyTitle) : ModeBase(dim,{rf},cf,keyTitle) {}

  virtual ~ModeBase() {}

};


/////
// 1
/////

/// Implements a free mode, that is, \f$H=-\delta\,a^\dagger a\f$ in a fully exact way, that is \f$\ket{\Psi(t)}=e^{i\delta\,t\,a^\dagger a}\ket{\Psi(0)}\f$ \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<typename AveragingType=mode::Averaged>
class Mode 
  : public mode::Exact, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  Mode(const mode::Pars&, AveragingConstructorParameters&&... );
};


/** \cond */

/// This is provided only for convenient use in maker functions. In this case unitary and “full” interaction picture coincide
template<typename AveragingType=mode::Averaged>
struct ModeUIP : Mode<AveragingType>
{
  template<typename... AveragingConstructorParameters>
  ModeUIP(const mode::Pars& p, AveragingConstructorParameters&&... a) : Mode<AveragingType>(p,a...) {}
};

/** \endcond */

/////
// 2
/////

/// Same as Mode, but without exact propagation \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<typename AveragingType=mode::Averaged>
class ModeSch
  : public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  ModeSch(const mode::Pars&, AveragingConstructorParameters&&... );
};


/////
// 3
/////

/// Implements a pumped mode, that is \f$H=-\delta\,a^\dagger a+i\lp\eta a^\dagger-\hermConj\rp\f$ in interaction picture defined by the first term \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<typename AveragingType=mode::Averaged>
class PumpedMode 
  : public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedMode(const mode::ParsPumped&, AveragingConstructorParameters&&... );
};


/** \cond */

template<typename AveragingType=mode::Averaged>
struct PumpedModeUIP : PumpedMode<AveragingType>
{
  template<typename... AveragingConstructorParameters>
  PumpedModeUIP(const mode::ParsPumped& p, AveragingConstructorParameters&&... a) : PumpedMode<AveragingType>(p,a...) {}
};

/** \endcond */

/////
// 4
/////

/// Same as PumpedMode, without exact propagation \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<typename AveragingType=mode::Averaged>
class PumpedModeSch
  : public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedModeSch(const mode::ParsPumped&, AveragingConstructorParameters&&... );
};


/////
// 5
/////

/// Implements a mode damped with rate \f$\kappa\f$, that is \f$H=\lp-\delta-i\kappa\rp a^\dagger a\f$, in a fully exact way, that is \f$\ket{\Psi(t)}=e^{-z\,t\,a^\dagger a}\ket{\Psi(0)}\f$, and \f$\Liou\rho=2\kappa\lp(n_\text{Th}+1)\,a\rho a^\dagger+n_\text{Th}\,a^\dagger\rho a\rp\f$ \see \ref genericelementsfreesmode "Summary of the various Mode classes"
/** \tparam TEMPERATURE governs whether the possibility of finite temperature is considered */
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class LossyMode 
  : public mode::Liouvillean<TEMPERATURE>, public mode::Exact, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  LossyMode(const mode::ParsLossy&, AveragingConstructorParameters&&... );

};

/////
// 6
/////

/// Same as LossyMode, but in unitary interaction picture, defined only by the \f$-\delta\,a^\dagger a\f$ part of the Hamiltonian \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class LossyModeUIP 
  : public mode::Liouvillean<TEMPERATURE>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  LossyModeUIP(const mode::ParsLossy&, AveragingConstructorParameters&&... );

};

/////
// 7
/////

/// Same as LossyMode, but in Schrödinger picture \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class LossyModeSch 
  : public mode::Liouvillean<TEMPERATURE>, public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  LossyModeSch(const mode::ParsLossy&, AveragingConstructorParameters&&... );

};


/////
// 8
/////

/// Combines LossyMode with pumping in full (non-unitary) interaction picture \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class PumpedLossyMode 
  : public mode::Liouvillean<TEMPERATURE>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyMode(const mode::ParsPumpedLossy&, AveragingConstructorParameters&&... );
};

/////
// 9
/////

/// Combines LossyModeUIP and PumpedMode \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class PumpedLossyModeUIP 
  : public mode::Liouvillean<TEMPERATURE>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyModeUIP(const mode::ParsPumpedLossy&, AveragingConstructorParameters&&... );
};

/////
// 10
/////

/// Combines LossyModeSch and PumpedModeSch \see \ref genericelementsfreesmode "Summary of the various Mode classes"
template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class PumpedLossyModeSch
  : public mode::Liouvillean<TEMPERATURE>, public mode::Hamiltonian<false>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyModeSch(const mode::ParsPumpedLossy&, AveragingConstructorParameters&&... );
};


//////////////////////////////////////////////////////////////////////
// One more to help testing the alternative jumping in MCWF_Trajectory
//////////////////////////////////////////////////////////////////////

template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class PumpedLossyModeAlternative 
  : public mode::Liouvillean<TEMPERATURE,true>, public mode::Hamiltonian<true>, public ModeBase, public AveragingType
{
public:
  template<typename... AveragingConstructorParameters>
  PumpedLossyModeAlternative(const mode::ParsPumpedLossy&, AveragingConstructorParameters&&... );
};


//////////////////////////////////////////////////////////////////////
// One more to test time-dependent jump and averages 
//////////////////////////////////////////////////////////////////////

class PumpedLossyModeIP_NoExact
  : public ModeBase,
    public quantumoperator::TridiagonalHamiltonian<1,true>, 
    public structure::ElementLiouvillean<1,1,true>,
    public structure::ElementAveraged<1,true>
{
public:
  PumpedLossyModeIP_NoExact(const mode::ParsPumpedLossy&);

private:
  using OneTime = structure::OneTime;

  void   doActWithJ (OneTime, mode::StateVectorLow&           ) const;
  double rate       (OneTime, const mode::LazyDensityOperator&) const;

  const mode::Averages average_v(OneTime, const mode::LazyDensityOperator&) const;

  const dcomp z_;

};

#endif // CPPQEDELEMENTS_FREES_MODE__H_INCLUDED
