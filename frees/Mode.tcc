// -*- C++ -*-
#ifndef ELEMENTS_FREES_IMPL_MODE_TCC_INCLUDED
#define ELEMENTS_FREES_IMPL_MODE_TCC_INCLUDED

#include "Mode_.h"

#include <boost/make_shared.hpp>


namespace mode {


#define TEMPLATE_PARAM_TEMP(temp) temp,A
#define SWITCH_helper(name,templateParam)				\
  using boost::make_shared;                                             \
  switch (qmp) {							\
  case QMP_IP  : return make_shared<name##Mode   <templateParam> >(p,a); \
  case QMP_UIP : return make_shared<name##ModeUIP<templateParam> >(p,a); \
  case QMP_SCH : ;							\
  }									\
  return make_shared<name##ModeSch<templateParam> >(p,a);


template<typename A>
const Ptr make(const Pars& p, QM_Picture qmp, const A& a)
{
  SWITCH_helper( ,A)
}


template<typename A>
const Ptr make(const ParsLossy& p, QM_Picture qmp, const A& a)
{
  if (!p.kappa) return make(static_cast<const Pars&>(p),qmp,a);
  else {
    if (p.nTh) { SWITCH_helper(Lossy,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(Lossy,TEMPLATE_PARAM_TEMP(false) ) }
  }
}


template<typename A>
const Ptr make(const ParsPumped& p, QM_Picture qmp, const A& a)
{
  if (!isNonZero(p.eta)) return make(static_cast<const Pars&>(p),qmp,a);
  else { SWITCH_helper(Pumped,A) }
}


template<typename A>
const Ptr make(const ParsPumpedLossy& p, QM_Picture qmp, const A& a)
{
  if      (!p.kappa)          return make(static_cast<const ParsPumped&>(p),qmp,a);
  else if (!isNonZero(p.eta)) return make(static_cast<const ParsLossy &>(p),qmp,a);
  else { 
    if (p.nTh) { SWITCH_helper(PumpedLossy,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(PumpedLossy,TEMPLATE_PARAM_TEMP(false) ) }
  }
}

#undef  SWITCH_helper
#undef  TEMPLATE_PARAM_TEMP


template<typename Base>
AveragedMonitorCutoff<Base>::AveragedMonitorCutoff()
  : Base(boost::assign::list_of("|Psi(cutoff-1)|^2"),KeyLabels())
{
}


template<typename Base>
const typename AveragedMonitorCutoff<Base>::Averages AveragedMonitorCutoff<Base>::average_v(const LazyDensityOperator& matrix) const
{
  const Averages averagesFromBase(Base::average_v(matrix));
  Averages averages(averagesFromBase.size()+1);
  averages(blitz::Range(0,averages.size()-2))=averagesFromBase;
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

#define BASE_init(R,C) ModeBase(p.cutoff,(R),(C))
#define BASE_initR(R)  BASE_init((R),ComplexFreqs())
#define BASE_initC(C)  BASE_init(RealFreqs(),(C))

#define TUPLE_delta(ISIP) "delta",-imag(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),BOOST_PP_IF(ISIP,1,p.cutoff)
#define TUPLE_eta "eta",get_eta(),sqrt(p.cutoff)
#define TUPLE_kappadelta(ISIP) "(kappa*(2*nTh+1),delta)",conj(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),BOOST_PP_IF(ISIP,1,p.cutoff)
#define TUPLE_kappa "kappa*(2*nTh+1)",real(get_zSch()),p.cutoff


template<typename A>
Mode<A>::Mode(const mode::Pars& p, const A& a) 
  : mode::Exact(dcomp(0,-p.delta),p.cutoff),
    BASE_initR(FREQS(TUPLE_delta(1))),
    A(a)
{}


template<typename A>
ModeSch<A>::ModeSch(const mode::Pars& p, const A& a) 
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),0,p.cutoff),
    BASE_initR(FREQS(TUPLE_delta(0))),
    A(a)
{
  getParsStream()<<"# Schroedinger picture.\n";
}


template<typename A>
PumpedMode<A>::PumpedMode(const mode::ParsPumped& p, const A& a)
  : mode::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta,p.cutoff),
    BASE_init(FREQS(TUPLE_delta(1)),FREQS(TUPLE_eta)),
    A(a)
{
  getParsStream()<<"# Pumped.\n";
}


template<typename A>
PumpedModeSch<A>::PumpedModeSch(const mode::ParsPumped& p, const A& a)
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),p.eta,p.cutoff),
    BASE_init(FREQS(TUPLE_delta(0)),FREQS(TUPLE_eta)),
    A(a)
{
  getParsStream()<<"# Pumped, Schroedinger picture.\n";
}


template<bool IS_FINITE_TEMP, typename A>
LossyMode<IS_FINITE_TEMP,A>::LossyMode(const mode::ParsLossy& p, const A& a)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
    mode::Exact(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.cutoff),
    BASE_initC(FREQS(TUPLE_kappadelta(1))),
    A(a)
{
  getParsStream()<<"# Lossy."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
LossyModeUIP<IS_FINITE_TEMP,A>::LossyModeUIP(const mode::ParsLossy& p, const A& a)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),0.),dcomp(0.,-p.delta),0,p.cutoff),
    BASE_initR(FREQS(TUPLE_kappa)(TUPLE_delta(1))),
    A(a)
{
  getParsStream()<<"# Lossy, Unitary interaction picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
LossyModeSch<IS_FINITE_TEMP,A>::LossyModeSch(const mode::ParsLossy& p, const A& a)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),0,p.cutoff),
    BASE_initR(FREQS(TUPLE_kappa)(TUPLE_delta(0))),
    A(a)
{
  getParsStream()<<"# Lossy, Schroedinger picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
PumpedLossyMode<IS_FINITE_TEMP,A>::PumpedLossyMode(const mode::ParsPumpedLossy& p, const A& a)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh), 
    mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
    BASE_initC(FREQS(TUPLE_kappadelta(1))(TUPLE_eta)),
    A(a)
{
  getParsStream()<<"# PumpedLossy."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
PumpedLossyModeUIP<IS_FINITE_TEMP,A>::PumpedLossyModeUIP(const mode::ParsPumpedLossy& p, const A& a)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh), 
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),0.),dcomp(0.,-p.delta),p.eta,p.cutoff),
    BASE_init(FREQS(TUPLE_kappa)(TUPLE_delta(1)),FREQS(TUPLE_eta)),
    A(a)
{
  getParsStream()<<"# PumpedLossy, Unitary interaction picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
PumpedLossyModeSch<IS_FINITE_TEMP,A>::PumpedLossyModeSch(const mode::ParsPumpedLossy& p, const A& a)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh), 
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
    BASE_initC(FREQS(TUPLE_kappadelta(0))(TUPLE_eta)),
    A(a)
{
  getParsStream()<<"# PumpedLossy, Schroedinger picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


//////////////////////////////////////////////


template<typename A>
PumpedLossyModeAlternative<A>::PumpedLossyModeAlternative(const mode::ParsPumpedLossy& p, const A& a)
  : mode::Liouvillean<false,true>(p.kappa,p.nTh), 
    mode::Hamiltonian<true>(0,dcomp(p.kappa,-p.delta),p.eta,p.cutoff),
    BASE_initC(FREQS(TUPLE_kappadelta(1))(TUPLE_eta)),
    A(a)
{
  getParsStream()<<"# PumpedLossy---Alternative jumping.\n";
}


#undef  TUPLE_kappa
#undef  TUPLE_kappadelta
#undef  TUPLE_eta
#undef  TUPLE_delta

#undef  BASE_initC
#undef  BASE_initR
#undef  BASE_init


#endif // ELEMENTS_FREES_IMPL_MODE_TCC_INCLUDED
