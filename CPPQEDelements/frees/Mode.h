// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"

#include "MultiDiagonal.h"

#include "Pars.h"

namespace mode {

using namespace ::structure; using namespace ::quantumdata;
using MultiDiagonal = ::quantumoperator::MultiDiagonal<1> ;


MultiDiagonal aOp(size_t cutoff);

MultiDiagonal aDagOp(size_t cutoff);

MultiDiagonal nOp(size_t cutoff);


double photonnumber(lazy_density_operator<1> auto matrix)
{
  double res=0;
  for (size_t n=0; n<matrix.extents[0]; ++n) res+=n*_(matrix,n);
  return res;
}


auto hamiltonian(size_t cutoff, dcomp z, double omegaKerr, dcomp eta)
{
  std::string label{"Mode monolithic MultiDiagonal:"};
  MultiDiagonal ham, a{aOp(cutoff)}, aDag{aDagOp(cutoff)};

  if (abs(z)) {label+=" linear term;"; ham-=z*nOp(cutoff);}
  if (omegaKerr) {label+=" Kerr term;"; ham+=(omegaKerr/1.i) * (aDag|aDag|a|a);}
  if (abs(eta)) {label+=" drive term;"; ham+=twoTimesImagPartOf(eta*aDag);}

  return makeHamiltonianElement<1>(label,std::move(ham));
    //-z*nOp(cutoff) + (omegaKerr/1.i) * (aDag|aDag|a|a) + twoTimesImagPartOf(eta*aDag)
}



TimeIndependentJump<1> aJump(double fact);
TimeIndependentJump<1> aDagJump(double fact);

TimeIndependentSuperoperator<1> aSuperoperator(double fact);
TimeIndependentSuperoperator<1> aDagSuperoperator(double fact);


Lindblad<1> photonLoss(double kappa, double nTh);
Lindblad<1> photonGain(double kappa, double nTh);


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho)
{
  double pn=0, pnn=0; dcomp a=0;
  for (size_t n=1; n<rho.extents[0]; ++n) {
    pn +=n*  _(rho,n);
    pnn+=n*n*_(rho,n);
    a+=sqrt(n)*_(rho,n,n-1);
  }
  return hana::make_tuple(pn,pnn,a);
};


constexpr auto postProcessor(decltype(expectationValues) t) {
  return [] (std::invoke_result_t<decltype(t),StateVectorConstView<1>> & t) {
    using namespace hana::literals; t[1_c] -= sqr(t[0_c]);
  };
}


constexpr ::cppqedutils::LogTree label(decltype(expectationValues)) { return {"photon number","photon number variance","ladder operator"}; }


static_assert(::structure::expectation_values_ns::labelled_and_with_nonlinear_postprocessing<decltype(expectationValues),1>);


auto make(size_t cutoff, double delta, double omegaKerr, dcomp eta, double kappa, double nTh)
{
  Liouvillian<1> liouvillian;
  std::vector<SystemFrequencyDescriptor> freqs;

  dcomp z{kappa*(2*nTh+1),-delta};

  if (delta) freqs.emplace_back("δ",delta,cutoff);

  if (omegaKerr) freqs.emplace_back("ω",omegaKerr,sqr(cutoff));

  if (abs(eta)) freqs.emplace_back("η",eta,sqrt(cutoff));

  if (kappa) {
    liouvillian.push_back(photonLoss(kappa,nTh));
    freqs.emplace_back("κ",kappa,(nTh+1)*cutoff);
    if (nTh) liouvillian.push_back(photonGain(kappa,nTh));
  }

  return QuantumSystemDynamics { std::move(freqs), std::move(liouvillian), hamiltonian(cutoff,z,omegaKerr,eta), exact_propagator_ns::noOp, expectationValues };
}



struct Pars
{
  size_t cutoff;
  double delta, omegaKerr, kappa, nTh;
  dcomp eta;

  Pars(popl::OptionParser& op, std::string mod="")
  {
    using ::parameters::_;
    add(mod,op,"Mode",
        _("cutoff","Fock space cutoff",10,cutoff),
        _("delta","detuning",-10.,delta),
        _("omegaKerr","Kerr constant",0.,omegaKerr),
        _("eta","drive amplitude",dcomp(0),eta),
        _("kappa","decay rate",10.,kappa),
        _("nTh","thermal photon number",0.,nTh));
  }

    // minitFock(p.add<size_t>("minitFock",mod,"Mode initial Fock state",0)),
    // minit(p.add<dcomp>("minit",mod,"Mode initial field",0)),

};


auto make(const Pars& p)
{
  return make(p.cutoff,p.delta,p.omegaKerr,p.eta,p.kappa,p.nTh);
}


/*


UnaryDiagonalPropagator<> propagator(dcomp z);


TimeIndependentSuperoperator<1> sigmaSuperoperator(double gamma_m);
TimeIndependentSuperoperator<1> sigmaPlusSuperoperator(double gamma_p);
TimeIndependentSuperoperator<1> sigma_zSuperoperator(double gamma_p);

TimeIndependentRate<1> sigmaRate(double gamma_m);
TimeIndependentRate<1> sigmaPlusRate(double gamma_p);

Lindblad<1> loss(double gamma_m) {return {"loss", sigmaJump(gamma_m), sigmaRate(gamma_m), sigmaSuperoperator(gamma_m)};}
Lindblad<1> gain(double gamma_p) {return {"gain", sigmaPlusJump(gamma_p), sigmaPlusRate(gamma_p), sigmaPlusSuperoperator(gamma_p)};}
// Lindblad<1> dephasing(double gamma_phi) {return {"dephasing", sigma_zJump(gamma_phi), sigma_zRate(gamma_phi), sigma_zSuperoperator(gamma_phi)};}


static constexpr auto expectationValues = [] (lazy_density_operator<1> auto rho) { return hana::make_tuple(real(_(rho,0,0)),_(rho,0,1)); };

::cppqedutils::LogTree label(decltype(expectationValues)) { return {"population","polarization"}; }


StateVector<1> state0();// {return mode::fock(0,2);}
StateVector<1> state1();// {return mode::fock(1,2);}
StateVector<1> init(dcomp psi1);
*/

} // mode
