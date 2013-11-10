#include "Qbit_.h"

#include "ParsQbit.h"
#include "StateVector.tcc"
#include "ElementLiouvillean.tcc"

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>


using namespace mathutils; using boost::make_shared;


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

const Averaged::Averages Averaged::average_v(NoTime, const LazyDensityOperator& matrix) const
{
  Averages averages(4);
  averages=matrix(0),matrix(1),real(matrix(1,0)),imag(matrix(1,0));
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
  psi(0)= fact*psi(0);
  psi(1)=-fact*psi(1);
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
  getParsStream()<<"# Schroedinger picture.\n";
}

   
PumpedQbit::PumpedQbit(const qbit::ParsPumped& p)
  : qbit::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta),
    QbitBase(TUPLE_delta(1),TUPLE_eta)
{
  getParsStream()<<"# Pumped.\n";
}


PumpedQbitSch::PumpedQbitSch(const qbit::ParsPumped& p)
  : qbit::Hamiltonian<false>(dcomp(0,-p.delta),p.eta),
    QbitBase(TUPLE_delta(0),TUPLE_eta)
{
  getParsStream()<<"# Pumped, Schroedinger picture.\n";
}


LossyQbit::LossyQbit(const qbit::ParsLossy& p)
  : qbit::Liouvillean(p.gamma),
    Exact(dcomp(p.gamma,-p.delta)),
    QbitBase{TUPLE_gammadelta(1)}
{
  getParsStream()<<"# Lossy.\n";
}


LossyQbitSch::LossyQbitSch(const qbit::ParsLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<false>(dcomp(p.gamma,-p.delta),0),
    QbitBase{TUPLE_gamma,TUPLE_delta(0)}
{
  getParsStream()<<"# Lossy, Schroedinger picture.\n";
}


LossyQbitUIP::LossyQbitUIP(const qbit::ParsLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),0),
    QbitBase{TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"# Lossy, Unitary interaction picture.\n";
}


PumpedLossyQbit::PumpedLossyQbit(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(0,dcomp(p.gamma,-p.delta),p.eta),
    QbitBase{TUPLE_gammadelta(1),TUPLE_eta}
{
  getParsStream()<<"# PumpedLossy.\n";
}


PumpedLossyQbitUIP::PumpedLossyQbitUIP(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),p.eta),
    QbitBase({TUPLE_gamma,TUPLE_delta(1)},TUPLE_eta)
{
  getParsStream()<<"# PumpedLossy, Unitary interaction picture.\n";
}


PumpedLossyQbitSch::PumpedLossyQbitSch(const qbit::ParsPumpedLossy& p)
  : qbit::Liouvillean(p.gamma),
    qbit::Hamiltonian<false>(dcomp(p.gamma,-p.delta),p.eta),
    QbitBase{TUPLE_gammadelta(0),TUPLE_eta}
{
  getParsStream()<<"# PumpedLossy, Schroedinger picture.\n";
}


LossyQbitWithPhaseNoise::LossyQbitWithPhaseNoise(const qbit::ParsLossy& p, double gamma_parallel)
  : Exact(dcomp(p.gamma,-p.delta)), // gamma_parallel does not contribute to the Hamiltonian
    qbit::LiouvilleanPhaseNoise(p.gamma,gamma_parallel),
    QbitBase(RF{"gamma_parallel",gamma_parallel,1},TUPLE_gammadelta(1))
{
  getParsStream()<<"# LossyWithPhaseNoise.\n";
}


LossyQbitWithPhaseNoiseUIP::LossyQbitWithPhaseNoiseUIP(const qbit::ParsLossy& p, double gamma_parallel)
  : qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),0),
    qbit::LiouvilleanPhaseNoise(p.gamma,gamma_parallel),
    QbitBase{RF{"gamma_parallel",gamma_parallel,1},TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"# LossyWithPhaseNoise, Unitary interaction picture.\n";
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
    return boost::make_shared<PumpedLossyQbit   >(p);
  case QMP_UIP :
    return boost::make_shared<PumpedLossyQbitUIP>(p);
  case QMP_SCH :
    ;
  }
  return make_shared<PumpedLossyQbitSch>(p);
}


const Tridiagonal sigmadagsigmaop()
{
  return mode::nop(make_shared<QbitBase>());
}


const StateVector init(const dcomp& psi1)
{
  StateVector res(2);
  res()(0)=sqrt(1-sqrAbs(psi1)); res()(1)=psi1;
  return res;
}


} // qbit
