// -*- C++ -*-
#ifndef _MODE_ELEMENT_IMPLEMENTATION_H
#define _MODE_ELEMENT_IMPLEMENTATION_H

#include<boost/assign/list_of.hpp>
#include<boost/assign/std/list.hpp>

namespace mode {


#define TEMPLATE_PARAM_TEMP(temp) temp,A
#define SWITCH_helper(name,templateParam)				\
  switch (qmp) {							\
  case QMP_IP  : return SmartPtr(new name##Mode<templateParam>(p));	\
  case QMP_UIP : return SmartPtr(new name##ModeUIP<templateParam>(p));	\
  case QMP_SCH : ;							\
  }									\
  return SmartPtr(new name##ModeSch<templateParam>(p));


template<typename A>
const SmartPtr maker(const Pars& p, QM_Picture qmp, const A&)
{
  SWITCH_helper( ,A)
}


template<typename A>
const SmartPtr maker(const ParsLossy& p, QM_Picture qmp, const A& a)
{
  if (!p.kappa) return maker(static_cast<const Pars&>(p),qmp,a);
  else {
    if (p.nTh) { SWITCH_helper(Lossy,TEMPLATE_PARAM_TEMP(true ) ) }
    else       { SWITCH_helper(Lossy,TEMPLATE_PARAM_TEMP(false) ) }
  }
}


template<typename A>
const SmartPtr maker(const ParsPumped& p, QM_Picture qmp, const A& a)
{
  if (!isNonZero(p.eta)) return maker(static_cast<const Pars&>(p),qmp,a);
  else { SWITCH_helper(Pumped,A) }
}


template<typename A>
const SmartPtr maker(const ParsPumpedLossy& p, QM_Picture qmp, const A& a)
{
  if      (!p.kappa)          return maker(static_cast<const ParsPumped&>(p),qmp,a);
  else if (!isNonZero(p.eta)) return maker(static_cast<const ParsLossy &>(p),qmp,a);
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
const typename AveragedMonitorCutoff<Base>::Averages AveragedMonitorCutoff<Base>::average(const LazyDensityOperator& matrix) const
{
  const Averages averagesFromBase(Base::average(matrix));
  Averages averages(averagesFromBase.size()+1);
  averages(blitz::Range(0,averages.size()-2))=averagesFromBase;
  averages(averages.size()-1)=matrix(matrix.getDimension()-1);
  return averages;
}


template<typename Base>
void AveragedMonitorCutoff<Base>::process(Averages& averages) const
{
  Averages ranged(averages(blitz::Range(0,averages.size()-2)));
  Averaged::process(ranged);
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
Mode<A>::Mode(const mode::Pars& p) 
  : mode::Exact(dcomp(0,-p.delta),p.cutoff),
    BASE_initR(boost::assign::tuple_list_of(TUPLE_delta(1))),
    A()
{}


template<typename A>
ModeSch<A>::ModeSch(const mode::Pars& p) 
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),0,p.cutoff),
    BASE_initR(boost::assign::tuple_list_of(TUPLE_delta(0))),
    A()
{
  getParsStream()<<"# Schroedinger picture.\n";
}


template<typename A>
PumpedMode<A>::PumpedMode(const mode::ParsPumped& p)
  : mode::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta,p.cutoff),
    BASE_init(boost::assign::tuple_list_of(TUPLE_delta(1)),boost::assign::tuple_list_of(TUPLE_eta)),
    A()
{
  getParsStream()<<"# Pumped.\n";
}


template<typename A>
PumpedModeSch<A>::PumpedModeSch(const mode::ParsPumped& p)
  : mode::Hamiltonian<false>(dcomp(0,-p.delta),p.eta,p.cutoff),
    BASE_init(boost::assign::tuple_list_of(TUPLE_delta(0)),boost::assign::tuple_list_of(TUPLE_eta)),
    A()
{
  getParsStream()<<"# Pumped, Schroedinger picture.\n";
}


template<bool IS_FINITE_TEMP, typename A>
LossyMode<IS_FINITE_TEMP,A>::LossyMode(const mode::ParsLossy& p)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
    mode::Exact(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.cutoff),
    BASE_initC(boost::assign::tuple_list_of(TUPLE_kappadelta(1))),
    A()
{
  getParsStream()<<"# Lossy."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
LossyModeUIP<IS_FINITE_TEMP,A>::LossyModeUIP(const mode::ParsLossy& p)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),0.),dcomp(0.,-p.delta),0,p.cutoff),
    BASE_initR(boost::assign::tuple_list_of(TUPLE_kappa)(TUPLE_delta(1))),
    A()
{
  getParsStream()<<"# Lossy, Unitary interaction picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
LossyModeSch<IS_FINITE_TEMP,A>::LossyModeSch(const mode::ParsLossy& p)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh),
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),0,p.cutoff),
    BASE_initR(boost::assign::tuple_list_of(TUPLE_kappa)(TUPLE_delta(0))),
    A()
{
  getParsStream()<<"# Lossy, Schroedinger picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
PumpedLossyMode<IS_FINITE_TEMP,A>::PumpedLossyMode(const mode::ParsPumpedLossy& p)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh), 
    mode::Hamiltonian<true>(0,dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
    BASE_initC(boost::assign::tuple_list_of(TUPLE_kappadelta(1))(TUPLE_eta)),
    A()
{
  getParsStream()<<"# PumpedLossy."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
PumpedLossyModeUIP<IS_FINITE_TEMP,A>::PumpedLossyModeUIP(const mode::ParsPumpedLossy& p)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh), 
    mode::Hamiltonian<true>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),0.),dcomp(0.,-p.delta),p.eta,p.cutoff),
    BASE_init(boost::assign::tuple_list_of(TUPLE_kappa)(TUPLE_delta(1)),boost::assign::tuple_list_of(TUPLE_eta)),
    A()
{
  getParsStream()<<"# PumpedLossy, Unitary interaction picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


template<bool IS_FINITE_TEMP, typename A>
PumpedLossyModeSch<IS_FINITE_TEMP,A>::PumpedLossyModeSch(const mode::ParsPumpedLossy& p)
  : mode::Liouvillean<IS_FINITE_TEMP>(p.kappa,p.nTh), 
    mode::Hamiltonian<false>(dcomp(mode::finiteTemperatureHamiltonianDecay(p,*this),-p.delta),p.eta,p.cutoff),
    BASE_initC(boost::assign::tuple_list_of(TUPLE_kappadelta(0))(TUPLE_eta)),
    A()
{
  getParsStream()<<"# PumpedLossy, Schroedinger picture."; mode::isFiniteTempStream(getParsStream(),p.nTh,*this);
}


//////////////////////////////////////////////


template<typename A>
PumpedLossyModeAlternative<A>::PumpedLossyModeAlternative(const mode::ParsPumpedLossy& p)
  : mode::Liouvillean<false,true>(p.kappa,p.nTh), 
    mode::Hamiltonian<true>(0,dcomp(p.kappa,-p.delta),p.eta,p.cutoff),
    BASE_initC(boost::assign::tuple_list_of(TUPLE_kappadelta(1))(TUPLE_eta)),
    A()
{
  getParsStream()<<"# PumpedLossy---Alternative jumping.\n";
}


#undef  TUPLE_kappa
#undef  TUPLE_kappadelta
#undef  TUPLE_eta
#undef  TUPLE_delta


#undef  BASE_initS
#undef  BASE_init


#endif // _MODE_ELEMENT_IMPLEMENTATION_H
