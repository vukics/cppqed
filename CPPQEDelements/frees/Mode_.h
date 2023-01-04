// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the \ref genericelementsfreesmode bundle (tackling the dynamics of a single harmonic-oscillator mode)}
#pragma once

#include "ParsMode.h"

#include "QM_Picture.h"

#include "Free.h"
#include "FreeExact.h"
#include "Liouvillian.h"
#include "TridiagonalHamiltonian.h"


class ModeBase;

/// Contains helpers for the \ref genericelementsfreesmode bundle
namespace mode {


using namespace structure::freesystem; using structure::NoTime;


struct diagonalFrequencies
{
  dcomp
    zODE,   // this goes into the ODE stepping
    zExact; // this is taken care of by exact propagation
  double omegaInt; // this is permanent interaction picture, must appear in averages
};


struct diagonalFrequenciesKerr : diagonalFrequencies
{
  double omegaKerr; // the coefficient of the Kerr term
};


double photonnumber(const LazyDensityOperator& matrix);


void aJump(StateVectorLow& psi, double fact, double omegaInt = 0.); // fact = sqrt(2.*kappa*(nTh+1))

void aDagJump(StateVectorLow& psi, double fact, double omegaInt = 0.); // fact = sqrt(2.*kappa*nTh)


void aSuperoperator(const DensityOperatorLow& rho, DensityOperatorLow& drhodt, double fact); // fact = 2.*kappa*(nTh+1)

void aDagSuperoperator(const DensityOperatorLow& rho, DensityOperatorLow& drhodt, double fact); // fact = 2.*kappa*nTh


structure::Lindblad<1> photonLoss(double kappaTimes_nThPlus1, double omegaInt = 0.); // no need for std::optional here, since there is a sensible default

structure::Lindblad<1> photonGain(double kappaTimes_nTh, double omegaInt = 0.);


structure::ElementLiouvillian<1> makeLiouvillian(double kappa, double nTh, double omegaInt = 0.);


structure::ExpectationValue<1> photonnumberEV_Variance;

structure::ExpectationValue<1> ladderOperatorEV;


// Quadrature variances are somewhat problematic since they involve all the expectation values


structure::ExpectationValue<1> monitorCutoff;


const structure::ElementExpectationValues<1> defaultExpectationValues{
  .label{"Mode"},
  .expectationValues{photonnumberEV_Variance,ladderOperatorEV}
};


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


const Ptr make(const Pars&, QM_Picture, const structure::ElementExpectationValues<1>& = defaultExpectationValues);

const Ptr make(const ParsDissipative&, QM_Picture, const structure::ElementExpectationValues<1>& = defaultExpectationValues);

const Ptr make(const ParsDriven&, QM_Picture, const structure::ElementExpectationValues<1>& = defaultExpectationValues);

const Ptr make(const ParsDrivenDissipative&, QM_Picture, const structure::ElementExpectationValues<1>& = defaultExpectationValues);



template<bool B>
inline std::ostream& isFiniteTempStream(std::ostream& os, double nTh)
{
  if constexpr (B) return os<<" Finite temperature.\nnTh="<<nTh<<std::endl;
  else return os<<std::endl;
}

template<bool B>
inline double finiteTemperatureHamiltonianDecay(const ParsDissipative& p)
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


} // mode



////////////////
//
// Highest level
//
////////////////


struct ModeBase : structure::Free
{
  ModeBase(size_t dim, const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs, const structure::ElementExpectationValues<1>& = mode::defaultExpectationValues);

  virtual ~ModeBase() {}
  
  
  std::optional<structure::ElementLiouvillian<1>> liouvillian;

  std::optional<structure::ElementExpectationValues<1>> expectationValues;

};


/////
// 1
/////

/// Implements a free mode, that is, \f$H=-\delta\,a^\dagger a\f$ in a fully exact way, that is \f$\ket{\Psi(t)}=e^{i\delta\,t\,a^\dagger a}\ket{\Psi(0)}\f$ \see \ref genericelementsfreesmode "Summary of the various Mode classes"
class Mode 
  : public mode::Exact, public ModeBase
{
public:
  template<typename... AveragingConstructorParameters>
  Mode(const mode::Pars&, AveragingConstructorParameters&&... );
};


/** \cond */

/// This is provided only for convenient use in maker functions. In this case unitary and “full” interaction picture coincide
struct ModeUIP : Mode
{
  template<typename... AveragingConstructorParameters>
  ModeUIP(const mode::Pars& p, AveragingConstructorParameters&&... a) : Mode(p,a...) {}
};

/** \endcond */

/////
// 2
/////

/// Same as Mode, but without exact propagation \see \ref genericelementsfreesmode "Summary of the various Mode classes"
class ModeSch
  : public mode::Hamiltonian<false>, public ModeBase
{
public:
  template<typename... AveragingConstructorParameters>
  ModeSch(const mode::Pars&, AveragingConstructorParameters&&... );
};


/////
// 3
/////

/// Implements a driven mode, that is \f$H=-\delta\,a^\dagger a+i\lp\eta a^\dagger-\hermConj\rp\f$ in interaction picture defined by the first term \see \ref genericelementsfreesmode "Summary of the various Mode classes"
class DrivenMode 
  : public mode::Hamiltonian<true>, public ModeBase
{
public:
  template<typename... AveragingConstructorParameters>
  DrivenMode(const mode::ParsDriven&, AveragingConstructorParameters&&... );
};


/** \cond */
struct DrivenModeUIP : DrivenMode
{
  template<typename... AveragingConstructorParameters>
  DrivenModeUIP(const mode::ParsDriven& p, AveragingConstructorParameters&&... a) : DrivenMode(p,a...) {}
};

/** \endcond */

/////
// 4
/////

/// Same as DrivenMode, without exact propagation \see \ref genericelementsfreesmode "Summary of the various Mode classes"
class DrivenModeSch
  : public mode::Hamiltonian<false>, public ModeBase
{
public:
  template<typename... AveragingConstructorParameters>
  DrivenModeSch(const mode::ParsDriven&, AveragingConstructorParameters&&... );
};



/*

//////////////////////////////////////////////////////////////////////
// One more to help testing the alternative jumping in MCWF_Trajectory
//////////////////////////////////////////////////////////////////////

template<bool TEMPERATURE=false, typename AveragingType=mode::Averaged>
class DrivenDissipativeModeAlternative 
  : public mode::Liouvillian<TEMPERATURE,true>, public mode::Hamiltonian<true>, public ModeBase
{
public:
  template<typename... AveragingConstructorParameters>
  DrivenDissipativeModeAlternative(const mode::ParsDrivenDissipative&, AveragingConstructorParameters&&... );
};


//////////////////////////////////////////////////////////////////////
// One more to test time-dependent jump and averages 
//////////////////////////////////////////////////////////////////////

class DrivenDissipativeModeIP_NoExact
  : public ModeBase,
    public quantumoperator::TridiagonalHamiltonian<1,true>, 
    public structure::ElementLiouvillian<1,1,true>,
    public structure::ElementAveraged<1,true>
{
public:
  DrivenDissipativeModeIP_NoExact(const mode::ParsDrivenDissipative&);

private:
  using OneTime = structure::OneTime;

  void   doActWithJ (OneTime, mode::StateVectorLow&           ) const;
  double rate       (OneTime, const mode::LazyDensityOperator&) const;

  const mode::Averages average_v(OneTime, const mode::LazyDensityOperator&) const;

  const dcomp z_;

};

*/
