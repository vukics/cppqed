// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDELEMENTS_FREES_MODE_TCC_INCLUDED
#define CPPQEDELEMENTS_FREES_MODE_TCC_INCLUDED

#include "Mode_.h"

#include "TridiagonalHamiltonian.h"


namespace mode {

#define TEMPLATE_PARAM_TEMP(temp) temp,AveragingType
#define SWITCH_helper(name,templateParam)                                 \
  switch (qmp) {                                                         \
    case QMP_IP  : return std::make_shared<name##Mode   <templateParam> >(p,a...); \
    case QMP_UIP : return std::make_shared<name##ModeUIP<templateParam> >(p,a...); \
  case QMP_SCH : ;                                                         \
  }                                                                         \
  return std::make_shared<name##ModeSch<templateParam> >(p,a...);


template<typename AveragingType, typename... AveragingConstructorParameters >
const Ptr make(const Pars& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  SWITCH_helper( ,AveragingType)
}


template<typename AveragingType, typename... AveragingConstructorParameters >
const Ptr make(const ParsDissipative& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if (!p.kappa) return make<AveragingType>(static_cast<const Pars&>(p),qmp,a...);
  else {
    if (p.nTh) { SWITCH_helper(Dissipative,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(Dissipative,TEMPLATE_PARAM_TEMP(false) ) }
  }
}


template<typename AveragingType, typename... AveragingConstructorParameters >
const Ptr make(const ParsDriven& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if (!isNonZero(p.eta)) return make<AveragingType>(static_cast<const Pars&>(p),qmp,a...);
  else { SWITCH_helper(Driven,AveragingType) }
}


template<typename AveragingType, typename... AveragingConstructorParameters >
const Ptr make(const ParsDrivenDissipative& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if      (!p.kappa)          return make<AveragingType>(static_cast<const ParsDriven&>(p),qmp,a...);
  else if (!isNonZero(p.eta)) return make<AveragingType>(static_cast<const ParsDissipative &>(p),qmp,a...);
  else { 
    if (p.nTh) { SWITCH_helper(DrivenDissipative,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(DrivenDissipative,TEMPLATE_PARAM_TEMP(false) ) }
  }
}

#undef  SWITCH_helper
#undef  TEMPLATE_PARAM_TEMP


template<typename Base>
AveragedMonitorCutoff<Base>::AveragedMonitorCutoff()
  : Base({{1,"|Psi(cutoff-1)|^2"}})
{
}


template<typename Base>
void AveragedMonitorCutoff<Base>::process_v(Averages& averages) const
{
  if (averages.size()>1) {
    Averages ranged(averages(blitz::Range(0,averages.size()-2)));
    Base::process_v(ranged);
  }
}



} // mode


////////////////
//
// Highest level
//
////////////////

#define TUPLE_delta(ISIP) RF{"delta",-imag(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),BOOST_PP_IF(ISIP,1,p.cutoff)}
#define TUPLE_eta CF{"eta",get_eta(),sqrt(p.cutoff)}
#define TUPLE_kappadelta(ISIP) CF{"(kappa*(2*nTh+1),delta)",conj(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),BOOST_PP_IF(ISIP,1,p.cutoff)}
#define TUPLE_kappa RF{"kappa*(2*nTh+1)",real(get_zSch()),p.cutoff}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
Mode<AveragingType>::Mode(const mode::Pars& p, AveragingConstructorParameters&&... a)
  : mode::Exact(dcomp(0,-p.delta),p.omegaKerr,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(1)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
ModeSch<AveragingType>::ModeSch(const mode::Pars& p, AveragingConstructorParameters&&... a)
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),0,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(0)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Schroedinger picture.\n";
}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
DrivenMode<AveragingType>::DrivenMode(const mode::ParsDriven& p, AveragingConstructorParameters&&... a)
  : mode::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(1),TUPLE_eta),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Driven.\n";
}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
DrivenModeSch<AveragingType>::DrivenModeSch(const mode::ParsDriven& p, AveragingConstructorParameters&&... a)
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(0),TUPLE_eta),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Driven, Schroedinger picture.\n";
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DissipativeMode<TEMPERATURE,AveragingType>::DissipativeMode(const mode::ParsDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE>(p.kappa,p.nTh),
    mode::Exact(dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),-p.delta),p.omegaKerr,p.cutoff),
    ModeBase(p.cutoff,TUPLE_kappadelta(1)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Dissipative."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DissipativeModeUIP<TEMPERATURE,AveragingType>::DissipativeModeUIP(const mode::ParsDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),0.),dcomp(0.,-p.delta),0,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappa,TUPLE_delta(1)}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Dissipative, Unitary interaction picture."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DissipativeModeSch<TEMPERATURE,AveragingType>::DissipativeModeSch(const mode::ParsDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),-p.delta),0,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappa,TUPLE_delta(0)}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Dissipative, Schroedinger picture."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DrivenDissipativeMode<TEMPERATURE,AveragingType>::DrivenDissipativeMode(const mode::ParsDrivenDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappadelta(1),TUPLE_eta}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"DrivenDissipative."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DrivenDissipativeModeUIP<TEMPERATURE,AveragingType>::DrivenDissipativeModeUIP(const mode::ParsDrivenDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),0.),dcomp(0.,-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappa,TUPLE_delta(1)},TUPLE_eta),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"DrivenDissipative, Unitary interaction picture."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DrivenDissipativeModeSch<TEMPERATURE,AveragingType>::DrivenDissipativeModeSch(const mode::ParsDrivenDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay<TEMPERATURE>(p),-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappadelta(0),TUPLE_eta}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"DrivenDissipative, Schroedinger picture."; mode::isFiniteTempStream<TEMPERATURE>(getParsStream(),p.nTh);
}


//////////////////////////////////////////////


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
DrivenDissipativeModeAlternative<TEMPERATURE,AveragingType>::DrivenDissipativeModeAlternative(const mode::ParsDrivenDissipative& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillian<TEMPERATURE,true>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(0,dcomp(p.kappa,-p.delta),p.eta,p.omegaKerr,p.omegaKerrAlter,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappadelta(1),TUPLE_eta}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"DrivenDissipative---Alternative jumping.\n";
}


#undef  TUPLE_kappa
#undef  TUPLE_kappadelta
#undef  TUPLE_eta
#undef  TUPLE_delta


#endif // CPPQEDELEMENTS_FREES_MODE_TCC_INCLUDED
