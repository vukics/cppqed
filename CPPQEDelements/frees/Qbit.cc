// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Qbit_.h"


using namespace cppqedutils; using namespace structure;


TimeIndependentTerm<1> qbit::diagonalH(dcomp z)
{
  return [=] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {dpsidt(1)-=z*psi(1);};
}

TimeIndependentTerm<1> qbit::offDiagonalH(dcomp eta)
{
  return [=] (StateVectorConstView<1> psi, StateVectorView<1> dpsidt) {dpsidt(0)+=-conj(eta)*psi(1); dpsidt(1)+=eta*psi(0);};
}

UnaryDiagonalPropagator<> qbit::propagator(dcomp z)
{
  return {2,[=] (double t, UnaryDiagonalPropagator<>::Diagonal& d) {d[0]=1; d[1]=exp(-z*t);}};
}


TimeIndependentJump<1> qbit::sigmaJump(double gamma_m) { return [=] (StateVectorView<1> psi) {psi(0)=sqrt(2.*gamma_m)*psi(1); psi(1)=0;}; }

TimeIndependentJump<1> qbit::sigmaPlusJump(double gamma_p) { return [=] (StateVectorView<1> psi) {psi(1)=sqrt(2.*gamma_p)*psi(0); psi(0)=0;}; }

TimeIndependentJump<1> qbit::sigma_zJump(double gamma_phi) { return [=] (StateVectorView<1> psi) {double fact=sqrt(2.*gamma_phi); psi(0)*=fact; psi(1)*=-fact;}; }


TimeIndependentSuperoperator<1> qbit::sigmaSuperoperator(double gamma_m)
{
  return [=] (DensityOperatorConstView<1> rho, DensityOperatorView<1> drhodt) {drhodt(0,0)+=2.*gamma_m*rho(1,1);};
}

TimeIndependentSuperoperator<1> qbit::sigmaPlusSuperoperator(double gamma_p)
{
  return [=] (DensityOperatorConstView<1> rho, DensityOperatorView<1> drhodt) {drhodt(1,1)+=2.*gamma_p*rho(0,0);};
}



/*

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



static constexpr auto dummyRate=[](const qbit::LazyDensityOperator&) {return -1;};

}


LiouvillianPhaseNoise::LiouvillianPhaseNoise(double gamma_perpendicular, double gamma_parallel) 
  : structure::ElementLiouvillianStrategies<1,2>({sigmaJump(gamma_perpendicular),sigma_zJump(gamma_parallel)},
                                                 {dummyRate,dummyRate},
                                                 "DissipativeQbitWithPhaseNoise",{"excitation loss","phase noise"})
{}

LiouvillianDrivenPhaseNoise::LiouvillianDrivenPhaseNoise(double gamma_perpendicular, double gamma_pump, double gamma_parallel)
  : structure::ElementLiouvillianStrategies<1,3>({sigmaJump(gamma_perpendicular),sigmaPlusJump(gamma_pump),sigma_zJump(gamma_parallel)},
                                                 {dummyRate,dummyRate,dummyRate},
                                                 "DissipativeQbitWithPhaseNoise",{"excitation loss","incoherent pump","phase noise"})
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

   
DrivenQbit::DrivenQbit(const qbit::ParsDriven& p)
  : qbit::Hamiltonian<true>(0,dcomp(0,-p.delta),p.eta),
    QbitBase(TUPLE_delta(1),TUPLE_eta)
{
  getParsStream()<<"Driven.\n";
}


DrivenQbitSch::DrivenQbitSch(const qbit::ParsDriven& p)
  : qbit::Hamiltonian<false>(dcomp(0,-p.delta),p.eta),
    QbitBase(TUPLE_delta(0),TUPLE_eta)
{
  getParsStream()<<"Driven, Schroedinger picture.\n";
}


DissipativeQbit::DissipativeQbit(double delta, double gamma)
  : qbit::Liouvillian(gamma),
    Exact(dcomp(gamma,-delta)),
    QbitBase{TUPLE_gammadelta(1)}
{
  getParsStream()<<"Dissipative.\n";
}


DissipativeQbitSch::DissipativeQbitSch(double delta, double gamma)
  : qbit::Liouvillian(gamma),
    qbit::Hamiltonian<false>(dcomp(gamma,-delta),0),
    QbitBase{TUPLE_gamma,TUPLE_delta(0)}
{
  getParsStream()<<"Dissipative, Schroedinger picture.\n";
}


DissipativeQbitUIP::DissipativeQbitUIP(double delta, double gamma)
  : qbit::Liouvillian(gamma),
    qbit::Hamiltonian<true>(dcomp(gamma,0),dcomp(0,-delta),0),
    QbitBase{TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"Dissipative, Unitary interaction picture.\n";
}


DrivenDissipativeQbit::DrivenDissipativeQbit(const qbit::ParsDrivenDissipative& p)
  : qbit::Liouvillian(p.gamma),
    qbit::Hamiltonian<true>(0,dcomp(p.gamma,-p.delta),p.eta),
    QbitBase{TUPLE_gammadelta(1),TUPLE_eta}
{
  getParsStream()<<"DrivenDissipative.\n";
}


DrivenDissipativeQbitUIP::DrivenDissipativeQbitUIP(const qbit::ParsDrivenDissipative& p)
  : qbit::Liouvillian(p.gamma),
    qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),p.eta),
    QbitBase({TUPLE_gamma,TUPLE_delta(1)},TUPLE_eta)
{
  getParsStream()<<"DrivenDissipative, Unitary interaction picture.\n";
}


DrivenDissipativeQbitSch::DrivenDissipativeQbitSch(const qbit::ParsDrivenDissipative& p)
  : qbit::Liouvillian(p.gamma),
    qbit::Hamiltonian<false>(dcomp(p.gamma,-p.delta),p.eta),
    QbitBase{TUPLE_gammadelta(0),TUPLE_eta}
{
  getParsStream()<<"DrivenDissipative, Schroedinger picture.\n";
}


DissipativeQbitWithPhaseNoise::DissipativeQbitWithPhaseNoise(double delta, double gamma, double gamma_parallel)
  : Exact(dcomp(gamma,-delta)), // gamma_parallel does not contribute to the Hamiltonian
    qbit::LiouvillianPhaseNoise(gamma,gamma_parallel),
    QbitBase(RF{"gamma_parallel",gamma_parallel,1},TUPLE_gammadelta(1))
{
  getParsStream()<<"DissipativeWithPhaseNoise.\n";
}


DissipativeQbitWithPhaseNoiseUIP::DissipativeQbitWithPhaseNoiseUIP(double delta, double gamma, double gamma_parallel)
  : qbit::Hamiltonian<true>(dcomp(gamma,0),dcomp(0,-delta),0),
    qbit::LiouvillianPhaseNoise(gamma,gamma_parallel),
    QbitBase{RF{"gamma_parallel",gamma_parallel,1},TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"DissipativeWithPhaseNoise, Unitary interaction picture.\n";
}

DrivenDissipativeQbitWithPhaseNoise::DrivenDissipativeQbitWithPhaseNoise(const qbit::ParsDrivenDissipativePhaseNoise& p)
  : qbit::Hamiltonian<true>(0,dcomp(p.gamma,-p.delta),p.eta),
    qbit::LiouvillianPhaseNoise(p.gamma,p.gamma_parallel),
    QbitBase{RF{"gamma_parallel",p.gamma_parallel,1},TUPLE_gammadelta(1),TUPLE_eta}
{
  getParsStream()<<"DrivenDissipativeWithPhaseNoise.\n";
}

DrivenDissipativeQbitWithPhaseNoiseUIP::DrivenDissipativeQbitWithPhaseNoiseUIP(const qbit::ParsDrivenDissipativePhaseNoise& p)
  : qbit::Hamiltonian<true>(dcomp(p.gamma,0),dcomp(0,-p.delta),p.eta),
    qbit::LiouvillianPhaseNoise(p.gamma,p.gamma_parallel),
    QbitBase{RF{"gamma_parallel",p.gamma_parallel,1},TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"DrivenDissipativeWithPhaseNoise, Unitary interaction picture.\n";
}


DissipativeQbitWithIncoherentPumpAndPhaseNoise::DissipativeQbitWithIncoherentPumpAndPhaseNoise(double delta, double gamma, double gamma_pump, double gamma_parallel)
  : Exact(dcomp(gamma,-delta)), // gamma_parallel does not contribute to the Hamiltonian
    qbit::LiouvillianDrivenPhaseNoise(gamma,gamma_pump,gamma_parallel),
    QbitBase(RF{"gamma_parallel",gamma_parallel,1},TUPLE_gammadelta(1))
{
  getParsStream()<<"DissipativeWithPhaseNoise.\n";
}


DissipativeQbitWithIncoherentPumpAndPhaseNoiseUIP::DissipativeQbitWithIncoherentPumpAndPhaseNoiseUIP(double delta, double gamma, double gamma_pump, double gamma_parallel)
  : qbit::Hamiltonian<true>(dcomp(gamma,0),dcomp(0,-delta),0),
    qbit::LiouvillianDrivenPhaseNoise(gamma,gamma_pump,gamma_parallel),
    QbitBase{RF{"gamma_parallel",gamma_parallel,1},TUPLE_gamma,TUPLE_delta(1)}
{
  getParsStream()<<"DissipativeWithPhaseNoise, Unitary interaction picture.\n";
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


Ptr make(const ParsDrivenDissipative& p, QM_Picture qmp)
{
  switch (qmp) {
  case QMP_IP  :
    if (p.gamma==0 && std::abs(p.eta)==0)
      return std::make_shared<Qbit            >(p);
    if (p.gamma==0)
      return std::make_shared<DrivenQbit      >(p);
    if (std::abs(p.eta)==0)
      return std::make_shared<DissipativeQbit       >(p);
    return std::make_shared<DrivenDissipativeQbit   >(p);
  case QMP_UIP :
    if (p.gamma==0 && std::abs(p.eta)==0)
      return std::make_shared<QbitUIP         >(p);
    if (p.gamma==0)
      return std::make_shared<DrivenQbitUIP   >(p);
    if (std::abs(p.eta)==0)
      return std::make_shared<DissipativeQbitUIP    >(p);
    return std::make_shared<DrivenDissipativeQbitUIP>(p);
  case QMP_SCH :
    ;
  }
  if (p.gamma==0 && std::abs(p.eta)==0)
    return std::make_shared<QbitSch         >(p);
  if (p.gamma==0)
    return std::make_shared<DrivenQbitSch   >(p);
  if (std::abs(p.eta)==0)
    return std::make_shared<DissipativeQbitSch    >(p);
  return std::make_shared<DrivenDissipativeQbitSch>(p);
}


Ptr make(const ParsDrivenDissipativePhaseNoise& p, QM_Picture qmp)
{
  if (p.gamma_parallel) {
    if (qmp==QMP_UIP) {
      if (std::abs(p.eta))
        return std::make_shared<DrivenDissipativeQbitWithPhaseNoiseUIP>(p);
      else
        return std::make_shared<      DissipativeQbitWithPhaseNoiseUIP>(p);
    }
    else {
      if (std::abs(p.eta))
        return std::make_shared<DrivenDissipativeQbitWithPhaseNoise>(p);
      else
        return std::make_shared<      DissipativeQbitWithPhaseNoise>(p);
    }
  }
  return make(static_cast<const ParsDrivenDissipative&>(p),qmp);
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
*/
