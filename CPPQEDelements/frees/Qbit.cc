// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Qbit_.h"

#include "ParsQbit.h"

#include <boost/bind.hpp>


using namespace cppqedutils; using std::make_shared;


namespace qbit {


///////////
//
// Averaged
//
///////////


Averaged::Averaged() 
  : Base(keyTitle,{"rho00","rho11","real(rho10=<sigma>)","imag(\")"})
{
}

const Averages Averaged::average_v(NoTime, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());
  averages=matrix(0),matrix(1),real(matrix(1)(0)),imag(matrix(1)(0));
  return averages;
}


namespace {


void sigmaJump(qbit::StateVectorLow& psi, double gamma_perpendicular)
{
  psi(0)=sqrt(2.*gamma_perpendicular)*psi(1);
  psi(1)=0;
}


void sigma_zJump(qbit::StateVectorLow& psi, double gamma_parallel)
{
  double fact=sqrt(2.*gamma_parallel);
  psi(0)*= fact; psi(1)*=-fact;
}


double dummyProba(const qbit::LazyDensityOperator&)
{
  return -1;
}


}


LiouvilleanPhaseNoise::LiouvilleanPhaseNoise(double gamma_perpendicular, double gamma_parallel) 
  : structure::ElementLiouvilleanStrategies<1,2>(JumpStrategies(bind(sigmaJump  ,_1,gamma_perpendicular),
                                                                bind(sigma_zJump,_1,gamma_parallel)),
                                                 JumpRateStrategies(dummyProba),
                                                 "LossyQbitWithPhaseNoise",{"excitation loss","phase noise"})
{}

} // qbit

////////////////
//
// Highest level
//
////////////////


QbitBase::QbitBase(const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs)
  : ModeBase(2,realFreqs,complexFreqs,qbit::keyTitle), Averaged()
{}

#define TUPLE_delta(ISIP) RF{"delta",-imag(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),1}
#define TUPLE_eta CF{"eta",get_eta(),1}
#define TUPLE_gammadelta(ISIP) CF{"(gamma,delta)",conj(BOOST_PP_IF(ISIP,get_zI,get_zSch)()),1}
#define TUPLE_gamma RF{"gamma",real(get_zSch()),1}

Qbit::Qbit(const qbit::Pars& p)
  : Exact(dcomp(0,-p.delta)),
    QbitBase{TUPLE_delta(1)}
{}


QbitSch::QbitSch(const qbit::Pars& p)
  : qbit::Hamiltonian<false>(dcomp(0,-p.delta),0),
    QbitBase{TUPLE_delta(0)}
{
  getParsStream()<<"Schroedinger picture.\n";
}

   
PumpedQbit::PumpedQbit(const qbit::ParsPumped& p)
  : qbit::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta),
    QbitBase(TUPLE_delta(1),TUPLE_eta)
{
  getParsStream()<<"Pumped.\n";
}


PumpedQbitSch::PumpedQbitSch(const qbit::ParsPumped& p)
  : qbit::Hamiltonian<false>(dcomp(0,-p.delta),p.eta),
    QbitBase(TUPLE_delta(0),TUPLE_eta)
{
  getParsStream()<<"Pumped, Schroedinger picture.\n";
}


LossyQbit::LossyQbit(double delta, double gamma)
  : qbit::Liouvillean(gamma),
    Exact(dcomp(gamma,-delta)),
    QbitBase{TUPLE_gammadelta(1)}
{
  getParsStream()<<"Lossy.\n";
}


LossyQbitSch::LossyQbitSch(double delta, double gamma)
  : qbit::Liouvillean(gamma),
    qbit::Hamiltonian<false>(dcomp(gamma,-delta),0),
    QbitBase{TUPLE_gamma,TUPLE_delta(0)}
{
  getParsStream()<<"Lossy, Schroedinger picture.\n";
}


LossyQbitUIP::LossyQbitUIP(double delta, double gamma)
  : qbit::Liouvillean(gamma),
    qbit::Hamiltonian<true>(dcomp(gamma,0),dcomp(0,-delta),0),
    QbitBase{TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"Lossy, Unitary interaction picture.\n";
}


PumpedLossyQbit::PumpedLossyQbit(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(0,dcomp(p.gamma,-p.delta),p.eta),
    QbitBase{TUPLE_gammadelta(1),TUPLE_eta}
{
  getParsStream()<<"PumpedLossy.\n";
}


PumpedLossyQbitUIP::PumpedLossyQbitUIP(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),p.eta),
    QbitBase({TUPLE_gamma,TUPLE_delta(1)},TUPLE_eta)
{
  getParsStream()<<"PumpedLossy, Unitary interaction picture.\n";
}


PumpedLossyQbitSch::PumpedLossyQbitSch(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<false>(dcomp(p.gamma,-p.delta),p.eta),
    QbitBase{TUPLE_gammadelta(0),TUPLE_eta}
{
  getParsStream()<<"PumpedLossy, Schroedinger picture.\n";
}


LossyQbitWithPhaseNoise::LossyQbitWithPhaseNoise(double delta, double gamma, double gamma_parallel)
  : Exact(dcomp(gamma,-delta)), // gamma_parallel does not contribute to the Hamiltonian
    qbit::LiouvilleanPhaseNoise(gamma,gamma_parallel),
    QbitBase(RF{"gamma_parallel",gamma_parallel,1},TUPLE_gammadelta(1))
{
  getParsStream()<<"LossyWithPhaseNoise.\n";
}


LossyQbitWithPhaseNoiseUIP::LossyQbitWithPhaseNoiseUIP(double delta, double gamma, double gamma_parallel)
  : qbit::Hamiltonian<true>(dcomp(gamma,0),dcomp(0,-delta),0),
    qbit::LiouvilleanPhaseNoise(gamma,gamma_parallel),
    QbitBase{RF{"gamma_parallel",gamma_parallel,1},TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"LossyWithPhaseNoise, Unitary interaction picture.\n";
}

PumpedLossyQbitWithPhaseNoise::PumpedLossyQbitWithPhaseNoise(const qbit::ParsPumpedLossyPhaseNoise& p)
  : qbit::Hamiltonian<true>(0,dcomp(p.gamma,-p.delta),p.eta),
    qbit::LiouvilleanPhaseNoise(p.gamma,p.gamma_parallel),
    QbitBase{RF{"gamma_parallel",p.gamma_parallel,1},TUPLE_gammadelta(1),TUPLE_eta}
{
  getParsStream()<<"PumpedLossyWithPhaseNoise.\n";
}

PumpedLossyQbitWithPhaseNoiseUIP::PumpedLossyQbitWithPhaseNoiseUIP(const qbit::ParsPumpedLossyPhaseNoise& p)
  : qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),p.eta),
    qbit::LiouvilleanPhaseNoise(p.gamma,p.gamma_parallel),
    QbitBase{RF{"gamma_parallel",p.gamma_parallel,1},TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"PumpedLossyWithPhaseNoise, Unitary interaction picture.\n";
}

#undef  TUPLE_gamma
#undef  TUPLE_gammadelta
#undef  TUPLE_eta
#undef  TUPLE_delta

//////////
//
// Helpers
//
//////////

namespace qbit {


Ptr make(const ParsPumpedLossy& p, QM_Picture qmp)
{
  switch (qmp) {
  case QMP_IP  :
    if (p.gamma==0 && std::abs(p.eta)==0)
      return std::make_shared<Qbit            >(p);
    if (p.gamma==0)
      return std::make_shared<PumpedQbit      >(p);
    if (std::abs(p.eta)==0)
      return std::make_shared<LossyQbit       >(p);
    return std::make_shared<PumpedLossyQbit   >(p);
  case QMP_UIP :
    if (p.gamma==0 && std::abs(p.eta)==0)
      return std::make_shared<QbitUIP         >(p);
    if (p.gamma==0)
      return std::make_shared<PumpedQbitUIP   >(p);
    if (std::abs(p.eta)==0)
      return std::make_shared<LossyQbitUIP    >(p);
    return std::make_shared<PumpedLossyQbitUIP>(p);
  case QMP_SCH :
    ;
  }
  if (p.gamma==0 && std::abs(p.eta)==0)
    return std::make_shared<QbitSch         >(p);
  if (p.gamma==0)
    return std::make_shared<PumpedQbitSch   >(p);
  if (std::abs(p.eta)==0)
    return std::make_shared<LossyQbitSch    >(p);
  return std::make_shared<PumpedLossyQbitSch>(p);
}


Ptr make(const ParsPumpedLossyPhaseNoise& p, QM_Picture qmp)
{
  if (p.gamma_parallel) {
    if (qmp==QMP_UIP) {
      if (std::abs(p.eta))
        return std::make_shared<PumpedLossyQbitWithPhaseNoiseUIP>(p);
      else
        return std::make_shared<      LossyQbitWithPhaseNoiseUIP>(p);
    }
    else {
      if (std::abs(p.eta))
        return std::make_shared<PumpedLossyQbitWithPhaseNoise>(p);
      else
        return std::make_shared<      LossyQbitWithPhaseNoise>(p);
    }
  }
  return make(static_cast<const ParsPumpedLossy&>(p),qmp);
}


const Tridiagonal sigmadagsigmaop()
{
  return mode::nop(make_shared<QbitBase>());
}


StateVector init(dcomp psi1)
{
  StateVector res(2);
  res(0)=sqrt(1-sqrAbs(psi1)); res(1)=psi1;
  return res;
}


} // qbit
