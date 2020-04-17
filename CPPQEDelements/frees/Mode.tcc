// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef CPPQEDELEMENTS_FREES_MODE_TCC_INCLUDED
#define CPPQEDELEMENTS_FREES_MODE_TCC_INCLUDED

#include "Mode_.h"

#include "TridiagonalHamiltonian.tcc"


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
const Ptr make(const ParsLossy& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if (!p.kappa) return make<AveragingType>(static_cast<const Pars&>(p),qmp,a...);
  else {
    if (p.nTh) { SWITCH_helper(Lossy,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(Lossy,TEMPLATE_PARAM_TEMP(false) ) }
  }
}


template<typename AveragingType, typename... AveragingConstructorParameters >
const Ptr make(const ParsPumped& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if (!isNonZero(p.eta)) return make<AveragingType>(static_cast<const Pars&>(p),qmp,a...);
  else { SWITCH_helper(Pumped,AveragingType) }
}


template<typename AveragingType, typename... AveragingConstructorParameters >
const Ptr make(const ParsPumpedLossy& p, QM_Picture qmp, AveragingConstructorParameters&&... a)
{
  if      (!p.kappa)          return make<AveragingType>(static_cast<const ParsPumped&>(p),qmp,a...);
  else if (!isNonZero(p.eta)) return make<AveragingType>(static_cast<const ParsLossy &>(p),qmp,a...);
  else { 
    if (p.nTh) { SWITCH_helper(PumpedLossy,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(PumpedLossy,TEMPLATE_PARAM_TEMP(false) ) }
  }
}

#undef  SWITCH_helper
#undef  TEMPLATE_PARAM_TEMP


template<typename Base>
AveragedMonitorCutoff<Base>::AveragedMonitorCutoff()
  : Base(KeyLabels(1,"|Psi(cutoff-1)|^2"))
{
}


template<typename Base>
const typename AveragedMonitorCutoff<Base>::Averages AveragedMonitorCutoff<Base>::average_v(NoTime t, const LazyDensityOperator& matrix) const
{
  auto averages{Base::average_v(t,matrix)}; // This is already of the correct size, since nAvr knows about the size updated by the derived class
  averages(averages.size()-1)=matrix(matrix.getDimension()-1);
  return averages;
}


template<typename Base>
void AveragedMonitorCutoff<Base>::process_v(Averages& averages) const
{
  Averages ranged(averages(blitz::Range(0,averages.size()-2)));
  Averaged::process_v(ranged);
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
  : mode::Exact(dcomp(0,-p.delta),p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(1)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
ModeSch<AveragingType>::ModeSch(const mode::Pars& p, AveragingConstructorParameters&&... a)
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),0,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(0)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Schroedinger picture.\n";
}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedMode<AveragingType>::PumpedMode(const mode::ParsPumped& p, AveragingConstructorParameters&&... a)
  : mode::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(1),TUPLE_eta),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Pumped.\n";
}


template<typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedModeSch<AveragingType>::PumpedModeSch(const mode::ParsPumped& p, AveragingConstructorParameters&&... a)
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),p.eta,p.cutoff),
    ModeBase(p.cutoff,TUPLE_delta(0),TUPLE_eta),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Pumped, Schroedinger picture.\n";
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
LossyMode<TEMPERATURE,AveragingType>::LossyMode(const mode::ParsLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
    mode::Exact(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.cutoff),
    ModeBase(p.cutoff,TUPLE_kappadelta(1)),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Lossy."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
LossyModeUIP<TEMPERATURE,AveragingType>::LossyModeUIP(const mode::ParsLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),0.),dcomp(0.,-p.delta),0,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappa,TUPLE_delta(1)}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Lossy, Unitary interaction picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
LossyModeSch<TEMPERATURE,AveragingType>::LossyModeSch(const mode::ParsLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),0,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappa,TUPLE_delta(0)}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"Lossy, Schroedinger picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedLossyMode<TEMPERATURE,AveragingType>::PumpedLossyMode(const mode::ParsPumpedLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappadelta(1),TUPLE_eta}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"PumpedLossy."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedLossyModeUIP<TEMPERATURE,AveragingType>::PumpedLossyModeUIP(const mode::ParsPumpedLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),0.),dcomp(0.,-p.delta),p.eta,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappa,TUPLE_delta(1)},TUPLE_eta),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"PumpedLossy, Unitary interaction picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedLossyModeSch<TEMPERATURE,AveragingType>::PumpedLossyModeSch(const mode::ParsPumpedLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE>(p.kappa,p.nTh),
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappadelta(0),TUPLE_eta}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"PumpedLossy, Schroedinger picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


//////////////////////////////////////////////


template<bool TEMPERATURE, typename AveragingType> template<typename... AveragingConstructorParameters>
PumpedLossyModeAlternative<TEMPERATURE,AveragingType>::PumpedLossyModeAlternative(const mode::ParsPumpedLossy& p, AveragingConstructorParameters&&... a)
  : mode::Liouvillean<TEMPERATURE,true>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(0,dcomp(p.kappa,-p.delta),p.eta,p.cutoff),
    ModeBase(p.cutoff,{TUPLE_kappadelta(1),TUPLE_eta}),
    AveragingType(std::forward<AveragingConstructorParameters>(a)...)
{
  getParsStream()<<"PumpedLossy---Alternative jumping.\n";
}


#undef  TUPLE_kappa
#undef  TUPLE_kappadelta
#undef  TUPLE_eta
#undef  TUPLE_delta


#endif // CPPQEDELEMENTS_FREES_MODE_TCC_INCLUDED
